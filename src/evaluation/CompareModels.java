package evaluation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import utility.CommonFunctions;
import utility.Helper;

/**
 * Compare models from 2 folders
 * @author Tuan
 *
 */
public class CompareModels {

	public static void main(String[] args) throws Exception{
		
		String folder1 = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr1_5kb/large_domains_output/";
		String folder2 = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr1_5kb/localModels/";
		String folderResult = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr1_5kb/comparison_2consecutive_domains.txt";
		
		//compareFoldersOfConsecutiveDomainModels();
		
		compare5kbVs1MB();
		
		//compareFoldersModels();
	}
	
	public static void compareFoldersModels() throws Exception{
		for(int i = 1; i <= 23; i++){
			String chr = "chr" + i;
			if (i == 23) chr = "chrX";
			String folder1 = "C:/Users/Tuan/git/GenomeMDS/GenomeMDS/output/hsc/BM/50kb_bm/" + chr + "/";
			String folder2 = "C:/Users/Tuan/git/GenomeMDS/GenomeMDS/output/hsc/FL/50kb_fl/" + chr + "/";
			String outputResult = "C:/Users/Tuan/git/GenomeMDS/GenomeMDS/output/hsc/" + chr + ".txt";
			
			double cor = compare2ModelFolders(folder1, folder2, outputResult);
			System.out.println(chr + ":\t" + cor);
		}
	}
	
	public static void compareFoldersOfConsecutiveDomainModels() throws Exception{
		for(int i = 1; i <= 23; i++){
			String chr = "chr" + i;
			if (i == 23) chr = "chrX";
			String folder1 = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/new/" + chr + "_5kb/large_domains_output/";
			String folder2 = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/new/" + chr + "_5kb/localModels/";
			String outputResult = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/new/" + chr + "_5kb/comparison_2consecutive_domains.txt";
			
			compare2ModelFolders(folder1, folder2, outputResult);
		}
	}
	
	public static void compare5kbVs1MB() throws Exception{
		
		String modelFile1 = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr2_5kb/chr2_1493321403413_1mb.pdb";
		String mappingFile1 = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/chr2_5kb/chr2_mapping_1mb.txt";
		
		String modelFile2 = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/1mb/chr2/chr2_1mb_gm12878_list_1493244082287.pdb";
		String mappingFile2 = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/1mb/chr2/chr2_1mb_gm12878_list_coordinate_mapping.txt";
		String outputResult = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/compare5kb_vs_1mb.txt";
		
		double cor = 0.0;
		//double cor = compare2Models(modelFile1, mappingFile1, modelFile2, mappingFile2);
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outputResult, true)));
		
		pw.println("Correlation between 1MB vs. 1kb models");
		
		for(int i = 23; i <= 23; i++){
			String chr = "chr" + i;
			if (i == 23) chr = "chrX";
			String folder = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/" + chr + "_5kb";
			File fo = new File(folder);
			
			if (!fo.exists()) continue;
			
			for(File f : fo.listFiles()){
				if (f.getName().endsWith("_1mb.pdb")){
					
					modelFile1 = f.getAbsolutePath();
					mappingFile1 = folder + "/" + chr + "_mapping_1mb.txt";
					
					modelFile2 = CommonFunctions.find_pdb_file("C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/1mb/" + chr);
					mappingFile2 = CommonFunctions.find_mapping_file("C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/1mb/" + chr);
					
					cor = compare2Models(modelFile1, mappingFile1, modelFile2, mappingFile2);
					
					pw.println(modelFile1 + " vs. " + modelFile2 + ": " + cor);
					System.out.println(modelFile1 + " vs. " + modelFile2 + ": " + cor);
					
					break;					
				}				
			}
		}
				
		pw.close();
		
		System.out.println(cor);		
		
		
	}
	
	
	//compare 2 folder of models
	public static double compare2ModelFolders(String folder1, String folder2, String outputResult) throws Exception{
		
		File file1 = new File(folder1);
		
		DescriptiveStatistics ds = new DescriptiveStatistics();
		
		PrintWriter pw = new PrintWriter(outputResult);
		pw.println("Comparing models of 2 consecutive domains extracted from the chromosome models and reconstructed individually");
		
		double avg = 0;
		int count = 0;
		
		for(File f : file1.listFiles()){
			if (f.isDirectory()){
				String modelFile1 = CommonFunctions.find_pdb_file(f.getAbsolutePath());
				String mappingFile1 = CommonFunctions.find_mapping_file(f.getAbsolutePath());
				
				int regionID = Integer.parseInt(f.getName().replaceAll("region_", ""));
				
				String modelFile2 = folder2 + "/region_" + regionID + ".pdb";
				String mappingFile2 = folder2 + "/region_" + regionID + "_mapping.txt";
				
				if (!CommonFunctions.isExist(modelFile2)) continue;
				
				double cor = compare2Models(modelFile1, mappingFile1, modelFile2, mappingFile2);
				avg += cor;
				count++;
				
				pw.println(modelFile1 + " vs. " + modelFile2 + ": " + cor);
				
				if (!Double.isNaN(cor)){				
					ds.addValue(cor);
				}else{
					System.out.println(modelFile1 + " vs. " + modelFile2 + ": " + cor);
				}
				//System.out.println(cor);
				
				if (cor < 0.5){
					System.out.println(cor + "\t" + regionID);
				}
				
			}
		}
		
		System.out.printf("Min: %.2f, max: %.2f, mean: %.2f, median: %.2f", ds.getMin(), ds.getMax(), ds.getMean(), ds.getPercentile(50));
		System.out.printf("\nThreshold of percentile: 2: %.2f, 4: %.2f, 6: %.2f, 8: %.2f, 10: %.2f", ds.getPercentile(2), ds.getPercentile(4), ds.getPercentile(6), ds.getPercentile(8), ds.getPercentile(10));
		
		pw.printf("Min: %.2f, max: %.2f, mean: %.2f, median: %.2f", ds.getMin(), ds.getMax(), ds.getMean(), ds.getPercentile(50));
		pw.printf("\nThreshold of percentile: 2: %.2f, 4: %.2f, 6: %.2f, 8: %.2f, 10: %.2f", ds.getPercentile(2), ds.getPercentile(4), ds.getPercentile(6), ds.getPercentile(8), ds.getPercentile(10));
		
		pw.close();
		
		return avg/count;
		
	}
	
	/**
	 * Compare 2 models using their absolute coordinates
	 * @param modelFile1
	 * @param mappingFile1
	 * @param modelFile2
	 * @param mappingFile2
	 * @return
	 * @throws Exception
	 */
	public static double compare2Models(String modelFile1, String mappingFile1, String modelFile2, String mappingFile2) throws Exception{
		
		Helper helper = Helper.getHelperInstance();
		
		HashSet<Integer> commonRegion;
		HashMap<Integer,Integer> mapping1 = helper.loadCoordinateMapping(mappingFile1);
		HashMap<Integer,Integer> mapping2 = helper.loadCoordinateMapping(mappingFile2);
		
		commonRegion = new HashSet<Integer>(mapping1.values());
		commonRegion.retainAll(mapping2.values());
		double[][] model1 = helper.loadPDBStructure(modelFile1);
		double[][] model2 = helper.loadPDBStructure(modelFile2);
		
		List<Double> lstX1 = new ArrayList<Double>();
		List<Double> lstY1 = new ArrayList<Double>();
		List<Double> lstZ1 = new ArrayList<Double>();
		
		List<Double> lstX2 = new ArrayList<Double>();
		List<Double> lstY2 = new ArrayList<Double>();
		List<Double> lstZ2 = new ArrayList<Double>();
		
		for(int i = 0; i < model1.length; i++){
			if (commonRegion.contains(mapping1.get(i))){
				lstX1.add(model1[i][0]);
				lstY1.add(model1[i][1]);
				lstZ1.add(model1[i][2]);				
			}
		}
		
		for(int i = 0; i < model2.length; i++){
			if (commonRegion.contains(mapping2.get(i))){
				lstX2.add(model2[i][0]);
				lstY2.add(model2[i][1]);
				lstZ2.add(model2[i][2]);				
			}
		}
		
		
		double[][] dist1 = new double[lstX1.size()][lstX1.size()];
		for(int i = 0; i < lstX1.size(); i++){
			for(int j = i + 1; j < lstX1.size(); j++){
				double d = Math.sqrt(helper.calEuclidianDist(lstX1.get(i), lstY1.get(i), lstZ1.get(i), lstX1.get(j), lstY1.get(j), lstZ1.get(j)));
				dist1[i][j] = d;
				dist1[j][i] = d;
			}
		}
		
		double[][] dist2 = new double[lstX2.size()][lstX2.size()];
		for(int i = 0; i < lstX2.size(); i++){
			for(int j = i + 1; j < lstX2.size(); j++){
				double d = Math.sqrt(helper.calEuclidianDist(lstX2.get(i), lstY2.get(i), lstZ2.get(i), lstX2.get(j), lstY2.get(j), lstZ2.get(j)));
				dist2[i][j] = d;
				dist2[j][i] = d;
			}
		}
		
		double cor = Evaluate.calCorrelation(dist1, dist2);
		
		if (Double.isNaN(cor)){
			System.out.println(cor);
		}
		
		return cor;				
	}

}
