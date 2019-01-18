package utility;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import full_modeling.HierarchicalModeling;
import topological_domain.IdentifyDomains;
import valueObject.RegionVO;

public class ExtractLocalModel {

	public static void main(String[] args) throws Exception {
		
		String globalModel = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr1_5kb/chr1_1493312276468.pdb";
		String globalMapping = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr1_5kb/chr1_coordinate_mapping.txt";
		String domainFile = "E:/GM12878/loop_arrow_domain/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_whole.txt";
		String outputFolder = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr1_5kb/localModels";
		
		int chrId = HierarchicalModeling.chr_id;
		int res = HierarchicalModeling.resolution;
		
		
		for(int i = 1; i <= 1; i++){
			String chr = "chr" + i;
			if (i == 23) chr = "chrX";
			chrId = i;
			
			String folder = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/" + chr + "_5kb/";
			File fo = new File(folder);
			for(File f : fo.listFiles()){
				if (f.getName().startsWith(chr) && f.getName().endsWith(".pdb")){
					globalModel = f.getAbsolutePath();
					globalMapping = folder + "/" + chr + "_coordinate_mapping.txt";
					outputFolder = folder + "/localModels";
					CommonFunctions.delete_file(outputFolder);
					CommonFunctions.make_folder(outputFolder);
					
					extractLocalModel(globalModel, globalMapping, domainFile, outputFolder, chrId, res);
					break;
				}
			}
			
		}
		
		//extractLocalModel(globalModel, globalMapping, domainFile, outputFolder, chrId, res);

	}
	
	public static void extractLocalModel(String globalModel, String globalMapping, String domainFile, String outputFolder, int chr, int res)
						throws Exception{
		
		Helper helper = Helper.getHelperInstance();
		
		IdentifyDomains identifyDomain = new IdentifyDomains(domainFile, chr, res);
		List<RegionVO> regions = identifyDomain.get_all_regions();
		
		HashMap<Integer, Integer> mappingID2Cor = helper.loadCoordinateMapping(globalMapping);
		HashMap<Integer, Integer> mappingCor2ID = new HashMap<Integer,Integer>();
		for(int i : mappingCor2ID.keySet()){
			mappingID2Cor.put(mappingCor2ID.get(i), i);
		}
		
		double[][] str = helper.loadPDBStructure(globalModel);
		
		for(int t = 0; t < regions.size() - 1; t++){
			
			int start = regions.get(t).getStart();
			int end = regions.get(t + 1).getEnd();
			
			//List<Double> lstX = new ArrayList<Double>();
			//List<Double> lstY = new ArrayList<Double>();
			//List<Double> lstZ = new ArrayList<Double>();
			List<Double> lst = new ArrayList<Double>();
			
			String mappingFile = "region_" + t + "_mapping.txt";
			String modelFile = "region_" + t + ".pdb";
			PrintWriter pwMapping = new PrintWriter(outputFolder + "/" + mappingFile);
			
			HashMap<Integer,Integer> idToChr = new HashMap<Integer, Integer>();
			
			int k = 0;
			for(int i = 0; i < str.length; i++){
				if (mappingID2Cor.get(i) >= start && mappingID2Cor.get(i) < end){
					
					lst.add(str[i][0]);
					lst.add(str[i][1]);
					lst.add(str[i][2]);
					
					if (mappingID2Cor.get(i) < regions.get(t).getEnd()){
						idToChr.put(k, 0);
					}else{
						idToChr.put(k, 1);
					}
					
					pwMapping.println(mappingID2Cor.get(i) + "\t" + k);
					k++;
				}
			}
			pwMapping.close();
			
			double[] localStr = new double[lst.size()];
			for(int i = 0; i < lst.size(); i++){
				localStr[i] = lst.get(i);
			}
			
			helper.writeStructure(outputFolder + "/" + modelFile, localStr, idToChr, "", true);
		}
		
	}
	
	

}
