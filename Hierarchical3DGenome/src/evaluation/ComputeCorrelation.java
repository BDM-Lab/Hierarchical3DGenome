package evaluation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import utility.CommonFunctions;
import utility.Helper;
import valueObject.Constraint;

public class ComputeCorrelation {

	public static void main(String[] args) throws Exception {
		
		String modelFile = "C:/Users/Tuan/git/GenomeMDS/GenomeMDS/output/preb/preb_norm_100kb_chr1_1494365133384.pdb";
		String mappingFile = "C:/Users/Tuan/git/GenomeMDS/GenomeMDS/output/preb/preb_norm_100kb_chr1_coordinate_mapping.txt";
		String dataFile = "C:/Users/Tuan/git/GenomeMDS/GenomeMDS/input/preb_norm_100kb_chr1.txt";
		String outputFile = "C:/Users/Tuan/git/Hierarchical3DGenome/Hierarchical3DGenome/output/correlations.txt";
		
		
		//double cor1 = spearmanCorrelation(dataFile, modelFile, mappingFile);
		//System.out.println(cor1);
		
		
		if (args.length > 0){		
			outputFile = args[0];
		}
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile, true)));
		for(int i = 1; i <= 23; i++){
			String chr = "chr" + i;
			if (i == 23) chr = "chrX";
			dataFile = "input/" + chr + "_1kb_gm12878_list.txt";
			String modelFolder = "output/" + chr + "_1kb";
			File folder = new File(modelFolder);
			
			if (!folder.exists()) continue;
			
			modelFile = null;
			for(File f : folder.listFiles()){
				if (f.getName().startsWith(chr) && f.getName().endsWith(".pdb") && !f.getName().contains("1mb")){
					modelFile = f.getAbsolutePath();
					break;
				}
			}
			
			mappingFile = modelFolder + "/" + chr + "_coordinate_mapping.txt";
			
			if (modelFile != null){
				double cor = spearmanCorrelation(dataFile, modelFile, mappingFile);
				pw.println(modelFile + "\t" + dataFile + ": " + cor);
				System.out.println(cor);
			}
		}
		
		pw.close();
	}
	
	public static double spearmanCorrelation(String dataFile, String modelFile, String mappingFile) throws Exception{
		Helper helper = Helper.getHelperInstance();
		
		ArrayList<Integer> lstPos = new ArrayList<Integer>();
		List<Constraint> lstCons = helper.readContactList(dataFile, lstPos, 0.0);
		
		Map<Integer,Integer> mapIDToPos = helper.loadCoordinateMapping(mappingFile);
		Map<Integer,Integer> mapPosToID = new HashMap<Integer, Integer>();
		for(int i : mapIDToPos.keySet()){
			mapPosToID.put(mapIDToPos.get(i), i);
		}
		
		double[][] str = helper.loadPDBStructure(modelFile);
		
		int x, y, i = 0;
		double[] ifs = new double[lstCons.size()];
		double[] dists = new double[lstCons.size()];
		
		
		for(Constraint con: lstCons){			
			if (mapPosToID.keySet().contains(con.getPos1()) && mapPosToID.keySet().contains(con.getPos2())){
				
				x = mapPosToID.get(con.getPos1());
				y = mapPosToID.get(con.getPos2());
				
				ifs[i] = con.getIF();
				dists[i] = Math.sqrt(helper.calEuclidianDist(str[x][0], str[x][1], str[x][2], str[y][0], str[y][1], str[y][2]));
				i++;
			}
		}
		
		return new SpearmansCorrelation().correlation(ifs, dists);
		
	}

}
