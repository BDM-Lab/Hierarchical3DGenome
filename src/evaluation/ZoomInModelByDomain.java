package evaluation;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import full_modeling.HierarchicalModeling;
import topological_domain.IdentifyDomains;
import utility.Helper;
import valueObject.RegionVO;

/**
 * Zoom in a model by averaging points belonging the same domain
 * @author Tuan
 *
 */
public class ZoomInModelByDomain {

	public static void main(String[] args) throws Exception {
		String input_domain_file = "E:/GM12878/loop_arrow_domain/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_whole.txt";
		
		
		String model_file = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr2_5kb/chr2_1493321403413.pdb";
		String mapping_file = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr2_5kb/chr2_coordinate_mapping.txt";
		
		String output_file = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr2_5kb/chr2_1493321403413_1mb.pdb";
		String output_mapping = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr2_5kb/chr2_mapping_1mb.txt";
		//
		
		//zoomInByBinSize(model_file, mapping_file, output_file, output_mapping, (int)1e6);
		
		for(int i = 23; i <= 23; i++){
			String chr = "chr" + i;
			if (i == 23) chr = "chrX";
			String folder = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/" + chr + "_5kb";
			File fo = new File(folder);
			
			if (!fo.exists()) continue;
			
			for (File f : fo.listFiles()){
				if (f.getName().startsWith(chr + "_") && f.getName().endsWith(".pdb")){
					model_file = f.getAbsolutePath();
					mapping_file = folder + "/" + chr + "_coordinate_mapping.txt";
					output_file = model_file.replace(".pdb", "_1mb.pdb");
					output_mapping = folder + "/" + chr + "_mapping_1mb.txt";
					
					zoomInByBinSize(model_file, mapping_file, output_file, output_mapping, (int)1e6);
					break;
				}
			}
		}
	
	}
	
	public static void zoomInByBinSize(String inputModel, String mappingFile, String outputFile, String outputMapping, int res) throws Exception{
		Helper helper = Helper.getHelperInstance();
		
		double[][] str = helper.loadPDBStructure(inputModel);
		HashMap<Integer, Integer> mapping = helper.loadCoordinateMapping(mappingFile);
		
		ArrayList<Double> lstX = new ArrayList<Double>();
		ArrayList<Double> lstY = new ArrayList<Double>();
		ArrayList<Double> lstZ = new ArrayList<Double>();
		
		PrintWriter pw = new PrintWriter(outputMapping);
		
		int count = 0;
		int k = 1;
		double x=0, y=0, z=0;
		pw.println(0 + "\t" + 0);
		for(int i = 0; i < str.length; i++){
			if (mapping.get(i) < res * k){
				x += str[i][0];
				y += str[i][1];
				z += str[i][2];
				count++;
			}else{
				
				i--;
				if (count > 0){
					x /= count;
					y /= count;
					z /= count;
					lstX.add(x);
					lstY.add(y);
					lstZ.add(z);
					
					x = 0;
					y = 0;
					z = 0;
					
					pw.println(k * res + "\t" + lstX.size());
					
				}	
				k++;
				count = 0;
			}
		}
		
		if (count > 0){
			x /= count;
			y /= count;
			z /= count;
			lstX.add(x);
			lstY.add(y);
			lstZ.add(z);
			pw.println(k * res + "\t" + lstX.size());
		}
		
		pw.close();
		
		int n = lstX.size();
		double[] a = new double[n * 3];
		for(int i = 0; i < n; i++){
			a[i * 3] = lstX.get(i);
			a[i * 3 + 1] = lstY.get(i);
			a[i * 3 + 2] = lstZ.get(i);
		}
		
		helper.writeStructure(outputFile, a, null, "", true);
	}

	
	public static void zoomModelByDomain(String input_domain_file, String model_file, String mapping_file, String output_file, int chrId) throws Exception{
		Helper helper = Helper.getHelperInstance();
		
		IdentifyDomains domain_identifer = new IdentifyDomains(input_domain_file, chrId, HierarchicalModeling.resolution);		
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
		start = regions.get(id).getStart();
		end = regions.get(id).getEnd();
		
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
				start = regions.get(id).getStart();
				end = regions.get(id).getEnd();
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
	
	
	
	
	
	

}
