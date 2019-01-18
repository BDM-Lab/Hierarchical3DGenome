package utility;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import topological_domain.IdentifyDomains;
import valueObject.Constants;
import valueObject.Constraint;



public class DistanceEstimate {

	public static void main(String[] args) {
		

	}
	
	
	public static double distanceBewteen2Domains(List<Constraint> lstCons, double convertFactor, 
			int id1, int id2, String domainModelFolder) throws Exception{

		
		Helper helper = Helper.getHelperInstance();
		
		String domainFolder1 = domainModelFolder + "/region_" + id1;
		String domainFolder2 = domainModelFolder + "/region_" + id2;
		
		
		////////////////////////
		String model_file1 = CommonFunctions.find_pdb_file(domainFolder1);
		
		if (model_file1 == null){
			System.err.println("No model file:" + domainFolder1);
			return Double.NaN;
		}
		
		String mapping_file1 = CommonFunctions.find_mapping_file(domainFolder1);
		
		double[][] model1 = helper.loadPDBStructure(model_file1);
		
		if (model1.length == 0){
			System.err.println("Model file without a structure:" + domainFolder1);
			return Double.NaN;		
		}
		
		HashMap<Integer,Integer> mapping = helper.loadCoordinateMapping(mapping_file1);
		HashMap<Integer,Integer> mapping1 = new HashMap<Integer,Integer>();
		//reverse mapping
		for(Entry<Integer,Integer> en : mapping.entrySet()){
			mapping1.put(en.getValue(), en.getKey());			
		}
		
		int cx1=0,cy1=0,cz1=0;
		for(int i = 0; i < model1.length; i++){
			cx1 += model1[i][0];
			cy1 += model1[i][1];
			cz1 += model1[i][2];
		}
		cx1 /= model1.length;
		cy1 /= model1.length;
		cz1 /= model1.length;
		
		/////////////////////////////
		String model_file2 = CommonFunctions.find_pdb_file(domainFolder2);
		
		if (model_file2 == null) return Double.NaN;
		
		String mapping_file2 = CommonFunctions.find_mapping_file(domainFolder2);
				
		double[][] model2 = helper.loadPDBStructure(model_file2);
		
		if (model2.length == 0) return Double.NaN;
		
		mapping = helper.loadCoordinateMapping(mapping_file2);
		
		
		HashMap<Integer,Integer> mapping2 = new HashMap<Integer,Integer>();
		//reverse mapping
		for(Entry<Integer,Integer> en : mapping.entrySet()){
			mapping2.put(en.getValue(), en.getKey());			
		}
		
		int cx2=0,cy2=0,cz2=0;
		for(int i = 0; i < model2.length; i++){
			cx2 += model2[i][0];
			cy2 += model2[i][1];
			cz2 += model2[i][2];
		}
		cx2 /= model2.length;
		cy2 /= model2.length;
		cz2 /= model2.length;
		
		/////////////////////////////////////////
		
		double dist1, dist2, minDist = Integer.MAX_VALUE;
		int pos1, pos2;
		
		for(Constraint con : lstCons){
			
				
			if (!mapping1.containsKey(con.getPos1())) continue;
			pos1 = mapping1.get(con.getPos1());		
			
			dist1 = Math.sqrt(helper.calEuclidianDist(cx1,cy1,cz1, model1[pos1][0], model1[pos1][1], model1[pos1][2]));
			
			
			if (!mapping2.containsKey(con.getPos2())) continue;
			pos2 = mapping2.get(con.getPos2());			
			
			dist2 = Math.sqrt(helper.calEuclidianDist(cx2,cy2,cz2, model2[pos2][0], model2[pos2][1], model2[pos2][2]));
			
			if (dist1 + dist2 + Constants.AVG_DIST/Math.pow(con.getIF(), convertFactor) < minDist)
				minDist = dist1 + dist2 + Constants.AVG_DIST/Math.pow(con.getIF(), convertFactor);
			
		}		
		
		return minDist;		
	}
	
	public static double distanceEstimateBetweenDomains(List<Constraint> lstCons, double convertFactor, String lowResGlobalModel, 
			String domainFile, String domainModelFolder,  int chrId, int res, PrintWriter... pw)throws Exception{
		
		Helper helper = Helper.getHelperInstance();
		
		double[][] globalModel = helper.loadPDBStructure(lowResGlobalModel);
		
		DescriptiveStatistics ds = new DescriptiveStatistics();
		IdentifyDomains identifyDomain = new IdentifyDomains(domainFile, chrId, res);
		System.out.println();
		for(int i = 0; i < globalModel.length - 1; i++){
			List<Constraint> lstConsDomain = identifyDomain.extractContactBetweenDomain(lstCons, i, i + 1);			
			double dist =  DistanceEstimate.distanceBewteen2Domains(lstConsDomain, convertFactor, i, i + 1, domainModelFolder);
			
			if (Double.isNaN(dist)) continue;
			
			double d = Math.sqrt(helper.calEuclidianDist(globalModel[i][0], globalModel[i][1], globalModel[i][2], 
					globalModel[i+1][0], globalModel[i + 1][1], globalModel[i + 1][2]));
			
			double ratio = dist/d;
			ds.addValue(ratio);
			if (pw.length > 0){
				pw[0].printf(" %.2f ", ratio);
			}
			System.out.printf(" %.2f ", ratio);
		}
		
		System.out.println("\nMedian:" + ds.getPercentile(50) + "\tMin:" + ds.getMin() + "\tMax:" + ds.getMax());
		if (pw.length > 0){
			pw[0].println("\nMedian:" + ds.getPercentile(50) + "\tMin:" + ds.getMin() + "\tMax:" + ds.getMax());
		}
		return ds.getPercentile(50);
	}
	
	public static double[][] distanceEstimateBetweenAllDomains(List<Constraint> lstCons, List<Integer> lstPos, double convertFactor, 
			String domainFile, String domainModelFolder, int chrId, int res)throws Exception{
		
		Helper helper = Helper.getHelperInstance();
		
		///
		//mapping coordinates to ID
		HashMap<Integer, Integer> coordinatesToID = new HashMap<Integer, Integer>();
		Collections.sort(lstPos);
		for(int i = 0; i < lstPos.size(); i++){
			coordinatesToID.put(lstPos.get(i), i);
		}
		//////////////////
		
		
		IdentifyDomains identifyDomain = new IdentifyDomains(domainFile, chrId, res);
		int n = identifyDomain.get_all_regions().size();
		
		double[][] dist = new double[n][n];
		for(int i = 0; i < n; i++){
			for(int j = i + 1; j < n; j++){
				dist[i][j] = Double.MAX_VALUE;
				dist[j][i] = Double.MAX_VALUE;
			}
		}
		
		//for each local model, mapping local ID to global ID, calculate the distance of that ID to its center
		HashMap<Integer, Double> distToCenter = new HashMap<Integer, Double>();
		File domainFolder = new File(domainModelFolder);
		//File[] localFolders = domainFolder.listFiles();
		for(File f : domainFolder.listFiles()){
			if (f.isDirectory()){
				String modelFile = CommonFunctions.find_pdb_file(f.getAbsolutePath());		
				if (modelFile == null) continue;
				
				String mappingFile = CommonFunctions.find_mapping_file(f.getAbsolutePath());
				
				HashMap<Integer,Integer> mapping = helper.loadCoordinateMapping(mappingFile);
				
				double[][] model = helper.loadPDBStructure(modelFile);
				double cx = 0, cy = 0, cz = 0;
				for(int i = 0; i < model.length; i++){
					cx += model[i][0];
					cy += model[i][1];
					cz += model[i][2];
				}
				cx /= model.length;
				cy /= model.length;
				cz /= model.length;
				
				for(int i = 0; i < model.length; i++){
					double d = Math.sqrt(helper.calEuclidianDist(cx, cy, cz, model[i][0], model[i][1], model[i][2]));
					if (coordinatesToID.containsKey(mapping.get(i))){
						distToCenter.put(mapping.get(i), d);
					}
				}				
			}
		}		
	
		
		for(Constraint con:lstCons){
			double d = Constants.AVG_DIST/(Math.pow(con.getIF(), convertFactor));
			
			if (distToCenter.containsKey(con.getPos1()) && distToCenter.containsKey(con.getPos2()) &&
					dist[con.getDomainID1()][con.getDomainID2()] > d + distToCenter.get(con.getPos1()) + distToCenter.get(con.getPos2())){
				
				dist[con.getDomainID1()][con.getDomainID2()] = d + distToCenter.get(con.getPos1()) + distToCenter.get(con.getPos2());
				dist[con.getDomainID2()][con.getDomainID1()] = dist[con.getDomainID1()][con.getDomainID2()];
			}			
		}
		
		return dist;
		
		
	}	
	

}
