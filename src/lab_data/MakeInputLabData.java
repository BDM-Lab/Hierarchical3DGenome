package lab_data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import valueObject.Contact;

/**
 * this class reads Hi-C data from Lieberman's data, and generate input file for Genome Reconstructor,...
 * @author tuan
 *
 */
public class MakeInputLabData {

	private static int baseSize = 1000000; //resolution
	private static String res = "1mb";
	private static int chr;
	
	private static String inputFile = "input/lab_data/contact_ALL";
		
	private static String outputFolder = "output/primary/";
	private static String outputAsList = "";
	private static String outputAsMatrix = "";
	
	private static int HEADER_CHR1 = 1;
	private static int HEADER_POS1 = 2;
	private static int HEADER_CHR2 = 3;
	private static int HEADER_POS2 = 4;
	
	public static void main(String[] args) throws Exception{
		//I make this comment
		for(int i = 7; i <= 7 ;i++){
			intraContact(i);
		}
		
	}
	
	/**
	 * make intra-chromosomal contact for chromosome chr
	 * @param chr
	 */
	private static void intraContact(int chr) throws Exception{
		
		outputAsList = "IFList_Chr_" + chr + "_" + res +"_pri.txt"; 
		outputAsMatrix = "IFMatrix_Chr_" + chr + "_" + res +"_pri.txt";
				
		int max = 0, n = 0;

		max = getMaxIndex(chr,inputFile); 
		if (max > n){
			n = max;
		}
		
		
		//length is n + 1 because the last index is n;
		n++;
		//contain raw count
		//n = 501;
		double[][] a = new double[n][n];
		
		readRawMatrix(a, chr, inputFile);
		

		//writeOutputAsList(standardNorm(a),outputFolder + outputAsList );
		//writeOutputAsMatrix(standardNorm(a),outputFolder + outputAsMatrix );
		
		//writeOutputAsMatrix(ICENorm(a),outputFolder + outputAsMatrix );
		writeOutputAsList(ICENorm(a),outputFolder + outputAsList );
		
	}
	
	/**
	 * write out contacts as list for my chromosome reconstructor
	 * @param a
	 * @param fileName
	 * @throws Exception
	 */
	public static void writeOutputAsList(double[][] a, String fileName) throws Exception{
		
		PrintWriter pw = new PrintWriter(fileName);
		for(int i = 0; i < a.length; i++){
			for(int j = i + 1; j < a.length; j++){
				if (!Double.isNaN(a[i][j]) && a[i][j] > 0){
					pw.printf("%d \t %d \t %.6f \n",i+1, j+1, a[i][j]);

				}
			}
		}
		pw.close();

	}
	
	public static void writeOutputAsList(List<Contact> lst, String fileName) throws Exception{

		PrintWriter pw = new PrintWriter(fileName);
		for(Contact ct : lst){
			pw.printf("%d \t %d \t %.6f \n",ct.getPos1(), ct.getPos2(), ct.getIF());
		}
		pw.close();

	}

	
	/**
	 * write out contacts as matrix for Shrec3D
	 * @param a
	 * @param fileName
	 * @throws Exception
	 */
	public static void writeOutputAsMatrix(double[][] a, String fileName) throws Exception{
		double[] s = new double[a.length];
		for(int i = 0; i < a.length; i++){
			s[i] = sumVector(a[i]);
			if (Math.abs(s[i] - 0.0) < 0.00001){
				System.out.println(i);
			}
		}
		
		PrintWriter pw = new PrintWriter(fileName);
		for(int i = 0; i < a.length; i++){
			if (s[i] > 0.000001){
				for(int j = 0; j < a.length; j++){
					if (s[j] > 0.000001){
						pw.printf("%20.12f ",a[i][j]);
					}
				}
				pw.println();
			}
		}
		pw.close();
	}
	
	/**
	 * normalize using Iterative Correction 
	 * 
	 * https://liorpachter.wordpress.com/2013/11/17/imakaev_explained/
	 * 
    Calculate S_i=\sum_j W^k_{ij} and let \overline{S}_i = \frac{1}{n}\sum_{i=1}^n S_i.
    Set \Delta B^k_i = \frac{S_i}{\overline{S}_i}.
    Set W^{k+1}_{ij} = \frac{W^k_{ij}}{\Delta B_i \Delta B_j}.
    Set B^{k+1}_i = B^k_i \cdot \Delta B^k_i.


	 * @param a
	 * @return
	 */
	public static double[][] ICENorm(double[][] a){
		
		double total = 0.0;
		double[][] w = new double[a.length][a.length];
		
		for(int i = 0; i < a.length; i++){
			for(int j = 0; j < a.length; j++){
				w[i][j] = a[i][j];
			}
		}
		
		double[] s = new double[a.length];
		double[] b = new double[a.length];
		double mean;
		int count;
		
		for(int t = 0; t < 30; t++){
			
			for(int i = 0; i < w.length; i++){
				s[i] = sumVector(w[i]);
			}
			
			count = 0;
			for(int i = 0; i < w.length; i++){
				if (!Double.isNaN(s[i]) && s[i] > 0) {
					count++;
				}					
			}
			/*
			for(int i = 0; i < s.length; i++){
				System.out.printf("(%d, %.5f) , ", i,s[i]);
			}
			System.out.println();
			*/
			
			if (count == 0){
				System.out.println();
			}
			mean = sumVector(s) / count;
			//normalize by mean
			for(int i = 0; i < w.length; i ++){
				b[i] = s[i] / mean;
			}
			
			for(int i = 0; i < w.length; i++){
				for(int j = i+1; j < w.length; j++){
					
					if (b[i] <= 0 || b[j] <= 0) continue;
					
					w[i][j] = w[i][j] / (b[i] * b[j]);
					
					w[j][i] = w[i][j];
				}
			}			
		}
		
		//normalize elements so that sum of a row = 1
		//
		for(int i = 0; i < w.length; i++){
			s[i] = sumVector(w[i]);
		}
		
		
		for(int i = 0; i < w.length; i++){
			for(int j = i + 1; j < w.length; j++){
				if (s[i] <= 0) continue;
				w[i][j] /= s[i];
				w[j][i] = w[i][j];
			}
		}		
		for(int i = 0; i < a.length; i++){
			for(int j = i + 1; j < a.length; j++){
				if (!Double.isNaN(w[i][j])){
					total += w[i][j];		
				}
			}
		}		
		//scale normalized IF
		for(int i = 0; i < w.length; i++){
			for(int j = 0; j < w.length; j++){
				w[i][j] *= total;
			}
		}
		//
		
		
		return w;
	}
	
	/*
	 * @param a
	 * @return
	 */
	public static void ICENormList(List<Contact> lst){
		
		double total = 0.0;
		
		//double[][] w = new double[a.length][a.length];
		
		/*
		List<Contact> lstW = new ArrayList<Contact>();
		
		for(Contact ct : lst){
			total += ct.getIF();			
		}
		total /= 2;
		*/

		
		//double[] s = new double[a.length];
		//double[] b = new double[a.length];
		HashMap<Integer,Double> mapS = new HashMap<Integer,Double>();
		HashMap<Integer,Double> mapB = new HashMap<Integer,Double>();
		double mean=0.0;
		
		for(int t = 0; t < 50; t++){
			/*
			for(int i = 0; i < w.length; i++){
				s[i] = sumVector(w[i]);
			}*/
			mapS.clear();
			mapB.clear();
			mean = 0;
			for(Contact ct : lst){
				if (!mapS.containsKey(ct.getPos1())){
					mapS.put(ct.getPos1(), 0.0);
				}
				if (!mapS.containsKey(ct.getPos2())){
					mapS.put(ct.getPos2(), 0.0);
				}
				
				if (ct.getPos1() != ct.getPos2()){
					mapS.put(ct.getPos1(), mapS.get(ct.getPos1()) + ct.getIF());
					mapS.put(ct.getPos2(), mapS.get(ct.getPos2()) + ct.getIF());
					mean += ct.getIF() * 2;
					
				}else{
					mapS.put(ct.getPos1(), mapS.get(ct.getPos1()) + ct.getIF());					
					mean += ct.getIF();
				}
			}
			/*
			for(int i = 0; i < mapS.size(); i++){
				if (mapS.containsKey(i)){
					System.out.printf("(%d, %.5f) , ",i, mapS.get(i));
				}
			}
			System.out.println();
			*/
			mean = mean / mapS.size();
			
			/*/normalize by mean
			for(int i = 0; i < w.length; i ++){
				b[i] = s[i] / mean;
			}*/
			for (Entry<Integer, Double> entry : mapS.entrySet()){
				mapB.put(entry.getKey(), entry.getValue() / mean);
			}
			
			/*
			for(int i = 0; i < w.length; i++){
				for(int j = i + 1; j < w.length; j++){					
					w[i][j] = w[i][j] / (b[i] * b[j]);					
					w[j][i] = w[i][j];
				}
			}*/
			
			for(Contact ct : lst){
				if (ct.getPos1() != ct.getPos2()){
					ct.setIF(ct.getIF() / (mapB.get(ct.getPos1()) * mapB.get(ct.getPos2()) ));
				}
			}			
			
		}
		
		/*/normalize elements so that sum of a row = 1		
		for(int i = 0; i < w.length; i++){
			s[i] = sumVector(w[i]);
		}*/
		mapS.clear();
		for(Contact ct : lst){
			if (!mapS.containsKey(ct.getPos1())){
				mapS.put(ct.getPos1(), 0.0);
			}
			if (!mapS.containsKey(ct.getPos2())){
				mapS.put(ct.getPos2(), 0.0);
			}
			
			if(ct.getPos1() != ct.getPos2()){
				mapS.put(ct.getPos1(), mapS.get(ct.getPos1()) + ct.getIF());
				mapS.put(ct.getPos2(), mapS.get(ct.getPos2()) + ct.getIF());
			}else{
				mapS.put(ct.getPos1(), mapS.get(ct.getPos1()) + ct.getIF());
			}
		}
		
		/*
		for(int i = 0; i < w.length; i++){
			for(int j = i + 1; j < w.length; j++){
				w[i][j] /= s[i];
				w[j][i] = w[i][j];
			}
		}*/
		for(Contact ct : lst){
			ct.setIF(ct.getIF() / mapS.get(ct.getPos1()));			
		}
		
		/*/scale normalized IF
		for(int i = 0; i < w.length; i++){
			for(int j = 0; j < w.length; j++){
				w[i][j] *= total;
			}
		}
		/*/
		for(Contact ct : lst){		
			total += ct.getIF();
		}		
		for(Contact ct : lst){
			ct.setIF(ct.getIF() * total);			
		}
		
		
		
		//return lst;
	}

	
	/**
	 * 
	 * @param a: raw count
	 * @return normalized count matrix
	 */
	public static double[][] standardNorm(double[][] a){
		
		double[][] norm = new double[a.length][a.length];
		double[] s1 = new double[a.length];
		double total = 0;
		
		for(int i = 0; i < a.length; i++){
			s1[i] = sumVector(a[i]);
		}
		
		total = sumVector(s1) / 2;
		
		for(int i = 0; i < a.length; i++){
			for(int j = i + 1; j < a.length; j++){
				
				norm[i][j] = a[i][j] * total / (s1[i] * s1[j]);
				norm[j][i] = norm[i][j];
				
			}
		}
		
		return norm;
		
	}
	
	/**
	 * sum all elements of array
	 * @param a
	 * @return
	 */
	public static double sumVector(double[] a){
		double sum = 0.0;
		for(double d : a){
			if (!Double.isNaN(d)){
				sum = sum + d;
			}
		}
		return sum;
	}
	
	
	/**
	 * get the largest index fragment of the chromosome
	 * @param chr
	 * @param fileName
	 * @return
	 */
	private static int getMaxIndex(int chr, String fileName) throws Exception{
		File file = new File(fileName);
		FileReader fr = null;
		BufferedReader br = null;
		String ln;
		String[] st;
		int chr1, chr2, pos1, pos2;
		int max = 0;
		try{
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			
			Pattern splitRegex = Pattern.compile("[|\\s]+");
			while((ln = br.readLine()) != null){
				st = splitRegex.split(ln);
				
				if (st[1].equalsIgnoreCase("x")){
					st[1] = "23";
				}else if (st[1].equalsIgnoreCase("y")){
					st[1] = "24";
				}
				
				if (st[3].equalsIgnoreCase("x")){
					st[3] = "23";
				}else if (st[3].equalsIgnoreCase("y")){
					st[3] = "24";
				}
				
				chr1 = Integer.parseInt(st[HEADER_CHR1]);
				pos1 = Integer.parseInt(st[HEADER_POS1]);
				chr2 = Integer.parseInt(st[HEADER_CHR2]);
				pos2 = Integer.parseInt(st[HEADER_POS2]);
				
				if (chr1 == chr && chr2 == chr){
					pos1 /= baseSize;
					pos2 /= baseSize;
					
					if (pos1 > max){
						max = pos1;
					}
					if (pos2 > max){
						max = pos2;
					}
				}
				
			}
			
		}catch(Exception ex){
			ex.printStackTrace();
			throw ex;
		}finally{
			if (br != null){
				br.close();
			}
			if (fr != null){
				fr.close();
			}
		}
		
		return max;

	}
	
	/**
	 * read observed intra-chromosomal contact matrix
	 * @param a: contact matrix returned, new count will be added up into this matrix
	 * @param chr: chromosome that contact matrix will be calculated
	 * @param fileName: input file
	 */
	private static void readRawMatrix(double[][] a, int chr,String fileName) throws Exception{
		
		File file = new File(fileName);
		FileReader fr = null;
		BufferedReader br = null;
		String ln;
		String[] st;
		int chr1, chr2, pos1, pos2;
		try{
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			
			Pattern splitRegex = Pattern.compile("[|\\s]+");
			while((ln = br.readLine()) != null){
				st = splitRegex.split(ln);
				
				if (st[1].equalsIgnoreCase("x")){
					st[1] = "23";
				}else if (st[1].equalsIgnoreCase("y")){
					st[1] = "24";
				}
				
				if (st[3].equalsIgnoreCase("x")){
					st[3] = "23";
				}else if (st[3].equalsIgnoreCase("y")){
					st[3] = "24";
				}
				
				chr1 = Integer.parseInt(st[HEADER_CHR1]);
				pos1 = Integer.parseInt(st[HEADER_POS1]);
				chr2 = Integer.parseInt(st[HEADER_CHR2]);
				pos2 = Integer.parseInt(st[HEADER_POS2]);
				
				if (chr1 == chr && chr2 == chr ){
					pos1 /= baseSize;
					pos2 /= baseSize;

					//if (pos1 < 500 && pos2 < 500) {
					if (pos1 != pos2){
						a[pos1][pos2] ++;
						a[pos2][pos1] ++;
					}else{						
						//a[pos1][pos2] ++;
					}
					//}
				}
				
			}
			
		}catch(Exception ex){
			ex.printStackTrace();
			throw ex;
		}finally{
			if (br != null){
				br.close();
			}
			if (fr != null){
				fr.close();
			}
		}
				
		
		
		
	}
}
