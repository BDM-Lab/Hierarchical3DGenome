package utility;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import topological_domain.IdentifyDomains;

/**
 * Zoom in a model by averaging points belonging the same domain
 * @author Tuan
 *
 */
public class ZoomInModelByDomain {

	public static void main(String[] args) throws Exception {
		String input_domain_file = "E:/GM12878/loop_arrow_domain/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_whole.txt";
		String mapping_file = "C:/Users/Tuan/workspace/GenomeMDS/output/hierarchical_modeling/first_few_domains_2119999_coordinate_mapping.txt";
		String model_file = "C:/Users/Tuan/workspace/GenomeMDS/output/hierarchical_modeling/first_few_domains_2119999_1488218788767.pdb";
		String output_file = "C:/Users/Tuan/workspace/GenomeMDS/output/hierarchical_modeling/first_few_domains_2119999_1488218788767_reduced.pdb";
		int resolution = Constants.RESOLUTION;
		
	
	}
	
	public static void zoom_model_by_domain(String input_domain_file, String model_file, String mapping_file, String output_file, int chrId) throws Exception{
		Helper helper = Helper.getHelperInstance();
		
		IdentifyDomains domain_identifer = new IdentifyDomains(input_domain_file, chrId, Constants.RESOLUTION);		
		List<RegionVO> regions = domain_identifer.get_all_regions();
		
		PrintWriter pw = new PrintWriter("domain.txt");
		for(RegionVO reg : regions){
			pw.println(reg.getStart() + "\t" + reg.getEnd());
		}
		pw.close();
		
		HashMap<Integer,Integer> map_id_cor =  helper.loadCoordinateMapping(mapping_file);
		
		List<Double> lstX = new ArrayList<Double>();
		List<Double> lstY = new ArrayList<Double>();
		List<Double> lstZ = new ArrayList<Double>();
		
		double[][] str = helper.loadPDBStructure(model_file);
		
		double x = 0, y = 0, z = 0;
		int count = 0, start, end, id = 0;
		
		while(map_id_cor.get(0) > regions.get(id).getStart()) id++;
		id--;	
		start = regions.get(id).start;
		end = regions.get(id).end;
		
		for(int i = 0; i < str.length; i++){
			if (map_id_cor.get(i) >= end){
				
				if (count > 0){
					
					System.out.println(start + "\t" + end + "\t:" + count);
					
					x /= count;
					y /= count;
					z /= count;
					
					lstX.add(x);
					lstY.add(y);
					lstZ.add(z);
				}
				id++;
				start = regions.get(id).start;
				end = regions.get(id).end;
				x = 0;
				y = 0;
				z = 0;
				count = 0;
							
				
			}else if (map_id_cor.get(i) < start){
				System.err.println("Error, check it!");
			}else{
				x += str[i][0];
				y += str[i][1];
				z += str[i][2];
				count++;
			}
		}
		
		if (count > 0){
			x /= count;
			y /= count;
			z /= count;
			
			lstX.add(x);
			lstY.add(y);
			lstZ.add(z);
		}
		
		double[][] new_str = new double[lstX.size()][3];
		for(int i = 0; i < lstX.size(); i++){
			new_str[i][0] = lstX.get(i);
			new_str[i][1] = lstY.get(i);
			new_str[i][2] = lstZ.get(i);
		}
		
		helper.writeStructure(output_file, helper.makeMatrixForPDBOutput(new_str), null, "", true);	
	}
	
	public static void extract_model(String model_file, String mapping_file, String domain_file, int[] range, String output_file, int chrId) throws Exception{
		Helper helper = Helper.getHelperInstance();
		
		IdentifyDomains domain_identifer = new IdentifyDomains(domain_file, chrId, Constants.RESOLUTION);		
		List<RegionVO> regions = domain_identifer.get_all_regions();
		
		HashMap<Integer, Integer> map_id_cor= helper.loadCoordinateMapping(mapping_file);
		
		int start_id = 0, end_id = 0;
		for(int i = 0; i < regions.size(); i++){			
			if (regions.get(i).getStart() == range[0]) start_id = i;
			if (regions.get(i).getEnd() == range[1]) end_id = i;
		}
		
		for(int i : map_id_cor.keySet()){
			if (map_id_cor.get(i) == start_id){
				start_id = i;
				break;
			}
		}
		
		for(int i : map_id_cor.keySet()){
			if (map_id_cor.get(i) == end_id) {
				end_id = i;
				break;
			}
		}
		
		System.out.println("Model:" + model_file);
		helper.writeStructure(output_file, helper.makeMatrixForPDBOutput(helper.extractChr(helper.loadPDBStructure(model_file), start_id, end_id)), null, "", true);
		
		//HashMap<Integer,Integer> map_id_cor =  helper.loadCoordinateMapping(mapping_file);
		
	}
	
	
	
	
	

}
