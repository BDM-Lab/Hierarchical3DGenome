package utility;

import static utility.CommonFunctions.lookup_region;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import topological_domain.IdentifyDomains;

public class ExtractWithinAndConsecutiveDomainContact {
	
	//private static int chr_id = 1;
	//private static int resolution = 1000;
	
	public static void main(String[] args) throws Exception {
		
		String input_observed_data_file = "E:/GM12878/1kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_1kb.RAWobserved";//chr1_1kb_gm12878_list_0k_20000k, chr1_1kb_gm12878_list
		
		//String input_data_file = "E:/GM12878/1kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_1kb_gm12878_list.txt";//chr1_1kb_gm12878_list_0k_20000k, chr1_1kb_gm12878_list
		String output_file = "E:/GM12878/1kb_resolution_intrachromosomal/chr1/MAPQGE30/winthin_and_consecutive_domain_contact.txt";
		String input_domain_file = "E:/GM12878/loop_arrow_domain/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_whole.txt";
		String model_file = "C:/Users/Tuan/workspace/GenomeMDS/output/hierarchical_modeling/global/chr1_1kb_gm12878_list_domain_norm_1486179128401.pdb";
		String model_id_mapping_file = "C:/Users/Tuan/workspace/GenomeMDS/output/hierarchical_modeling/global/chr1_1kb_gm12878_list_domain_norm_coordinate_mapping.txt";
		
		//extract_within_and_consecutive_contact(input_data_file, input_domain_file, output_file);
		
		determine_proximal_domain_and_filter(input_observed_data_file, input_domain_file, model_file, model_id_mapping_file, output_file, 1);
				
	}
	
	public static void determine_proximal_domain_and_filter(String input_data_file, String input_domain_file, String model_file, 
			String model_id_mapping_file, String output_file, int chrId) throws Exception{
		Helper helper = Helper.getHelperInstance();
		double[][] str = helper.loadPDBStructure(model_file);
		
		HashMap<Integer, Integer> map_id_pos = helper.loadCoordinateMapping(model_id_mapping_file);
		
		//maximum genomic location id
		int max_pos = Collections.max(map_id_pos.values());
		
		//DescriptiveStatistics ds = new DescriptiveStatistics();
		ArrayList<Double> lst = new ArrayList<Double>();
		int n = str.length;
		double d = 0;
		for(int i = 0; i < str.length; i++){
			for(int j = i + 1; j < str.length; j++){
				d = Math.sqrt(helper.calEuclidianDist(str[i][0], str[i][1], str[i][2], str[j][0], str[j][1], str[j][2]));
				lst.add(d);
			}
		}
		
		
		boolean[][] is_proximal_domain = new boolean[max_pos + 1][max_pos + 1];
		for(int i = 0; i < max_pos + 1; i++){
			Arrays.fill(is_proximal_domain[i], false);
		}
		
		Collections.sort(lst);
		double median = lst.get(lst.size()/2);
		
		for(int i = 0; i < n - 1; i++){
			is_proximal_domain[i][i] = true;
			is_proximal_domain[i][i + 1] = true;
			is_proximal_domain[i + 1][i] = true;
		}
		is_proximal_domain[n - 1][n - 1] = true;
		
		/*
		for(int i = 0; i < n; i++){
			for(int j = i + 2; j < n; j++){//if |i - j| < 5, contacts between them are included
				//d = Math.sqrt(helper.calEuclidianDist(str[i][0], str[i][1], str[i][2], str[j][0], str[j][1], str[j][2]));
				//if (d > median){ // if two domains are too far away from each other
					is_proximal_domain[map_id_pos.get(i)][map_id_pos.get(j)] = false;
					is_proximal_domain[map_id_pos.get(j)][map_id_pos.get(i)] = false;
				//}
			}
		}
		*/
		
		extract_within_and_consecutive_contact(input_data_file, input_domain_file, output_file, chrId, is_proximal_domain);
	}
	
	public static void extract_within_and_consecutive_contact(String input_data_file, String input_domain_file, String output_file, int chrId, boolean[][] ... is_proximal_domain) throws Exception{
		
		IdentifyDomains domain_identifer = new IdentifyDomains(input_domain_file, chrId, Constants.RESOLUTION);
		
		List<RegionVO> regions = domain_identifer.get_all_regions();
		Collections.sort(regions);
		
		
		BufferedReader br = new BufferedReader(new FileReader(new File(input_data_file)));
		String ln;
		String[] st;
		RegionVO tmp = new RegionVO(1,0,0), domain1, domain2;
		
		int x, y;
		
		PrintWriter pw = new PrintWriter(output_file);
		while((ln = br.readLine()) != null) {
			st = ln.split("\\s+");
			x = Integer.parseInt(st[0]);
			y = Integer.parseInt(st[1]);
			//f = Double.parseDouble(st[2]);
			
			tmp.setChr_id(chrId);
			
			//look for domain that contains first index
			tmp.setStart(x);
			tmp.setEnd(x);			
			domain1= lookup_region(regions, tmp);
			
			//look for domain that contains second index
			tmp.setStart(y);
			tmp.setEnd(y);
			domain2 = lookup_region(regions, tmp);
			
			if (is_proximal_domain.length == 0){
				if (Math.abs(domain1.getId() - domain2.getId()) <= 1){
					pw.println(ln);
				}			
			}else{
				if (is_proximal_domain[0][domain1.getId()][domain2.getId()]){
					pw.println(ln);
				}
			}
		}
		pw.close();
		br.close();
	}

	
}
