package topological_domain;

import static utility.CommonFunctions.lookup_region;
import static utility.CommonFunctions.overlap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import utility.Constants;
import utility.RegionVO;
import valueObject.Constraint;


/**
 * This class currently works for one chromosome only, need upgrade for whole genome later
 * 
 * This class is to identify TAD from a Hi-C matrix input
 * @author Tuan
 *
 */
public class IdentifyDomains {
	
	String input_file;
	int resolution = Constants.RESOLUTION;
	int largest_sequence_id;
	int chrom;
	private List<RegionVO> domains = null;
	
	public static void main(String[] args)throws Exception{
		String input_file = "E:/GM12878/loop_arrow_domain/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_whole.txt";
		IdentifyDomains domain_identifer = new IdentifyDomains(input_file, 1, 1000);
		
		List<RegionVO> lst = domain_identifer.get_all_regions();	
		for(RegionVO reg : lst){
			System.out.println(reg.getChr_id() + "\t" + reg.getStart() + "\t" + reg.getEnd());
		}
	}
	
	
	
	/**
	 * 
	 * @param input_file: contact file
	 * @param res: resolution
	 * @throws Exception 
	 */
	public IdentifyDomains(String input_file, int chr, int res) throws Exception{
		this.input_file = input_file;
		this.resolution = res;
		this.chrom = chr;
		domains = get_all_regions();
			
	}
	
	
	public List<RegionVO> getDomains() {
		return domains;
	}

	public void setDomains(List<RegionVO> domains) {
		this.domains = domains;
	}
	
	//calculate coverage of all regions
	public long get_size(List<RegionVO> lst){
		Collections.sort(lst);
		
		long total = 0;
		int id = lst.get(0).getStart();
		int end = lst.get(lst.size() - 1).getEnd();
		int i = 0;
		while(id < end){
			if (lst.get(i).getStart() <= id && id < lst.get(i).getEnd()){
				total++;
				id++;
			}else if (id >= lst.get(i).getEnd()){
				i++;
			}else if (id < lst.get(i).getStart()){
				id = lst.get(i).getStart();
			}else{
				System.err.println("error in count coverage of regions");
			}			
		}
		
		return total;
	}
	
	public List<RegionVO> get_all_regions() throws Exception{
		determine_largest_sequence_id();
		
		List<RegionVO> rs = identify_domains();
		
		//long size1 = get_size(rs);
		
		rs = merge_domains(rs);
			
		//long size2 = get_size(rs);
		
//		if (size1 != size2){
//			System.err.println("error in merging domain");
//		}
//		
//		for(int i = 0; i < rs.size(); i++){
//			for(int j = i + 1; j < rs.size(); j++){
//				if (overlap(rs.get(i), rs.get(j)) > 0){
//					System.err.println("error in merging domain, they are still overlap");
//				}
//			}
//		}		
		
		List<RegionVO> regions =  add_non_domains(rs);
		
		regions = merge_short_domains(regions);
		
		//long size3 = get_size(regions);		
		//System.out.println("Coverage of the whole chromosome: " + size3);
		
		
//		//debug
//		HashMap<Integer, String> regSet = new HashMap<Integer, String>();
//		for(RegionVO reg : regions){
//			if (!regSet.containsKey(reg.getId())){				
//				regSet.put(reg.getId(), reg.getStart() + "-" + reg.getEnd());
//			}else{
//				System.out.println("Regions with same ID:" + reg.getId() + "\t" + regSet.get(reg.getId()) + "\t" + reg.getStart() + "-" + reg.getEnd());
//			}
//		}
//		//
		
		Collections.sort(regions);
		
		for(int i = 0; i < regions.size(); i++) regions.get(i).setId(i);
		
		return regions;
	}
	
	/**
	 * Length of chromosomes in the input data
	 */
	public void determine_largest_sequence_id(){
		this.largest_sequence_id = 249250612;//249,250,621;
	}
	
	public List<RegionVO> identify_domains() throws Exception{
		List<RegionVO> result = new ArrayList<RegionVO>();
		
		BufferedReader br = new BufferedReader(new FileReader(new File(input_file)));
		String ln;
		String[] st;
		int chr_id, x, y;
		RegionVO reg;
		while((ln = br.readLine()) != null){
			st = ln.split("\\s+");
			if (st[0].equals("chrX")) chr_id = 23;
			else chr_id = Integer.parseInt(st[0].replace("chr", ""));
			
			if (chr_id != chrom) continue; // testing with chromosome 1
			
			x = Integer.parseInt(st[1]);
			y = Integer.parseInt(st[2]);
			reg = new RegionVO(chr_id, x, y);
			reg.setDomain(true);
			result.add(reg);
		}
		
		br.close();	
		
		return result;
	}
	
	/*
	 * Merger domains that overlap or contain
	 */
	public List<RegionVO> merge_domains(List<RegionVO> rs){
		Collections.sort(rs);
		int i = 0, j = 1;
		RegionVO r1, r2, new_region;
		
		ArrayList<RegionVO> merged_list = new ArrayList<RegionVO>();
		merged_list.add(rs.get(0));
		
		while(j < rs.size()){
			r1 = merged_list.get(merged_list.size() - 1);
			r2 = rs.get(j);
			
			if (overlap(r1, r2) >= 1){
				new_region = new RegionVO(r1.getChr_id(), Math.min(r1.getStart(), r2.getStart()), Math.max(r1.getEnd(), r2.getEnd()));
				new_region.setDomain(true);
				//rs.set(i, new_region );		
				merged_list.set(merged_list.size() - 1, new_region);
			}else{
				//i++;
				//rs.set(i, r2);
				merged_list.add(r2);
			}			
			j++;			
		}
		
//		while(rs.size() > i + 1){
//			rs.remove(rs.size() - 1);
//		}
		
		return merged_list;
	}

	/**
	 * Merge too small domain( non-domains) into adjacent domains
	 * @param rs
	 * @return
	 */
	public List<RegionVO> merge_short_domains(List<RegionVO> rs){
		Collections.sort(rs);
		int i = 0, j = 1;
		RegionVO r1, r2, new_region;
		
		ArrayList<RegionVO> merged_list = new ArrayList<RegionVO>();
		merged_list.add(rs.get(0));
		
		while(j < rs.size()){
			r1 = merged_list.get(merged_list.size() - 1);
			r2 = rs.get(j);
			
			if ( ((r1.getEnd() - r1.getStart())/resolution < 300 && (r2.getEnd() - r2.getStart()) / resolution < 20 ) || (r1.getEnd() - r1.getStart())/resolution < 20){
				
				new_region = new RegionVO(r1.getChr_id(), Math.min(r1.getStart(), r2.getStart()), Math.max(r1.getEnd(), r2.getEnd()));
				new_region.setDomain(r1.isDomain() || r2.isDomain());

				merged_list.set(merged_list.size() - 1, new_region);
			}else{
				merged_list.add(r2);
			}			
			j++;			
		}	
		return merged_list;
	}

	
	/**
	 * Add regions between domains as non-domains, this works for one chromosome
	 * @param rs
	 */
	public List<RegionVO> add_non_domains(List<RegionVO> rs){
		List<RegionVO> new_list = new ArrayList<RegionVO>();
		RegionVO r1, r2, new_region;
		if (rs.size() > 0 && rs.get(0).getStart() > resolution){
			new_region = new RegionVO(rs.get(0).getChr_id(), 0, rs.get(0).getStart());
			new_list.add(new_region);
		}
		for(int i = 0; i < rs.size() - 1; i++){			
			r1 = rs.get(i);
			r2 = rs.get(i + 1);
			
			new_list.add(rs.get(i));
			
			//if there is a gap between r1 and r2, make a non-domain region between r1, r2
			if (r1.getChr_id() == r2.getChr_id() && r1.getEnd() + resolution < r2.getStart()){
				new_region = new RegionVO(r1.getChr_id(), r1.getEnd() + 1, r2.getStart() - 1);
				new_region.setDomain(false);
				new_list.add(new_region);
			}			
		}
		
		new_list.add(rs.get(rs.size() - 1));		
		
		r1 = rs.get(rs.size() - 1);
		if (r1.getEnd() + resolution < largest_sequence_id){
			new_region = new RegionVO(r1.getChr_id(), r1.getEnd() + 1, largest_sequence_id);
			new_list.add(new_region);
		}
		
		return new_list;
	}
	
	/**
	 * Assign domain for constraints
	 * @param lstCon
	 */
	public void assignDomainID(List<Constraint> lstCons, List<RegionVO> domains){
		RegionVO tmp = new RegionVO(1,0,0), domain1, domain2;
		int id1, id2, x, y;		
		
		int chr = domains.get(0).getChr_id();
		for(Constraint con : lstCons){
			
			x = con.getPos1();
			y = con.getPos2();
						
			tmp.setChr_id(chr);
			
			
				//look for domain that contains first index
			tmp.setStart(x);
			tmp.setEnd(x);			
			domain1= lookup_region(domains, tmp);
			id1 = domain1.getId();
			con.setDomainID1(id1);	
			
			tmp.setStart(y);
			tmp.setEnd(y);
			domain2 = lookup_region(domains, tmp);
			id2 = domain2.getId();
			con.setDomainID2(id2);
			
		}
	}
	
	//extract contact between domain1 and domain2
	public List<Constraint> extractContactBetweenDomain(List<Constraint> lstCons, int id1, int id2){
		ArrayList<Constraint> rs = new ArrayList<Constraint>();
		for(Constraint con:lstCons){
			if ((con.getDomainID1() == id1 && con.getDomainID2() == id2)
				|| (con.getDomainID1() == id2 && con.getDomainID2() == id1)) {
				rs.add(con);
			}
		}
		return rs;
	}
	
}
