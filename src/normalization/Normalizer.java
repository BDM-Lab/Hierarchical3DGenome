package normalization;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import lab_data.MakeInputLabData;
import valueObject.Constants;

/**
 * This is to package as an independent program
 * Perform normalization on the input data
 * @author Tuan
 *
 */
public class Normalizer {

	public static void main(String[] args) throws Exception {
		
		if (args.length < 6) {
			System.err.println("Input format error!");
			System.err.println("-i {input file} -o {ouput file} -m {method: 1: coverage normalization, 2: iterative correction of Hi-C}");
			System.exit(1);
		}
		String fileIn = null;
		String fileOut = null;
		int method = 1;
		for(int i = 0; i < args.length; i++){
			if (args[i].equals("-i")){
				fileIn = args[i+1];
			}else if (args[i].equals("-o")){
				fileOut = args[i+1];
			}else if (args[i].equals("-m")){
				method = Integer.parseInt(args[i+1].trim());
			}
		}
		
		if (fileIn == null || fileOut == null){
			System.err.println("Input format error!");
			System.err.println("-i {input file} -o {ouput file} -m {method: 1: coverage normalization, 2: iterative correction of Hi-C}");
		}
		
		normalize(fileIn, fileOut, method);
		
//		ArrayList<Integer> lstPos = new ArrayList<Integer>();
//		HashMap<Integer,Integer> hmPos = new HashMap<Integer,Integer>();
//		
//		System.out.println("Reading input ...");
//		double[][] a = readContactData(fileIn, lstPos,hmPos);
//		
//		System.out.println("Performing normalization ...");
//		if (method == 1){
//			a = MakeInputLabData.standardNorm(a);
//		}else{
//			a = MakeInputLabData.ICENorm(a);
//		}
//		
//		System.out.println("Writing out ...");
//		writeOut(fileOut,a,lstPos);
		
	}
	
	public static void normalize(String input_file, String output_file, int method) throws Exception{
		ArrayList<Integer> lstPos = new ArrayList<Integer>();
		HashMap<Integer,Integer> hmPos = new HashMap<Integer,Integer>();
		
		System.out.println("Reading input ...");
		double[][] a = readContactData(input_file, lstPos,hmPos);
		
		System.out.println("Performing normalization ...");
		if (method == 1){
			a = MakeInputLabData.standardNorm(a);
		}else{
			a = MakeInputLabData.ICENorm(a);
		}
		
		System.out.println("Writing out ...");
		writeOut(output_file,a,lstPos);
	}
	
	public static void normalize(String input_folder, int method) throws Exception{
		File folder = new File(input_folder);
		File[] files = folder.listFiles();
		
		int numOfcores = Runtime.getRuntime().availableProcessors();		
		numOfcores = Math.min(numOfcores * 2 , Constants.MAX_NUM_THREAD);
		
		ExecutorService executor = Executors.newFixedThreadPool(numOfcores);
		
		for(File f:files){
			if(f.length() > 0){			
				//normalize(f.getAbsolutePath(), f.getAbsolutePath(), method);
				executor.submit(new WorkerThread(f.getAbsolutePath(), f.getAbsolutePath()));
			}
		}		
		executor.shutdown();
		
		try{
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
		}catch(InterruptedException e){
			e.printStackTrace();
		}
	}
	static class WorkerThread implements Runnable {
		
		String input_file, output_file;
		public WorkerThread(String in_file, String out_file){
			this.input_file = in_file;
			this.output_file = out_file;
		}
		@Override
		public void run(){
			try {
				normalize(input_file, output_file, 1);
			} catch (Exception e) {
				e.printStackTrace();
			}			
		}		
	}
		
	
	public static void writeOut(String fileOut, double[][] a, ArrayList<Integer> lstPos) throws Exception{
		PrintWriter pw = new PrintWriter(fileOut);
		for(int i = 0; i < a.length; i++){
			for(int j = i + 1; j < a.length; j++){
				if (a[i][j] > 0){
					pw.println(lstPos.get(i) + "\t" + lstPos.get(j) + "\t" + a[i][j]);
				}
			}			
		}
		pw.flush();
		pw.close();
		
	}

	
	public static double[][] readContactData(String inputFile, ArrayList<Integer> lstPos, HashMap<Integer,Integer> hmPos )	throws FileNotFoundException, Exception{
		
		//contact matrix will be returned
		double[][] a;
		
		int x,y,id1,id2;
		double f;
		int nbr = 3; // number of numbers in each line
		
		File file = new File(inputFile);

		FileReader fr=null;
		BufferedReader br = null;		
		Pattern splitRegex = Pattern.compile("[:\\s]+");
		
		HashSet<Integer> posSet = new HashSet<Integer>();
		StringBuilder sb = new StringBuilder();
		try{
				
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			String ln;
			String[] st;
			int count = 1;
			ln = br.readLine();
			sb.append(ln).append(" ");
			nbr = splitRegex.split(ln.trim()).length;
			long progress = 0;
			while((ln = br.readLine()) != null){
				if (ln.trim().length() == 0 || !Character.isDigit(ln.charAt(0)) ){
					continue;
				}
				
				//read every 10 thousand lines and split at once
				sb.append(ln).append(" ");
				count++;
				
				if (count == 200000){
					count = 0;
					st = splitRegex.split(sb.toString());
					sb = new StringBuilder();
					
					if (st.length % nbr != 0){
						throw new Exception("There is a line that doesn't contain exactly 3 numbers");
					}
					//each line contains 'nbr' numbers
					for(int i = 0; i < st.length / nbr; i++){
						
						//only take first 3 numbers
						x = Integer.parseInt(st[i * nbr + 0]);						
						//position2
						y = Integer.parseInt(st[i * nbr + 1]);
						
						//interaction frequency
						//f = Double.parseDouble(st[i * nbr + 2]);
						
						//keeping absolute positions, so that later they can be recovered from indices
						
						posSet.add(x);
						
						posSet.add(y);
						 
					}
				}
				progress++;
				//System.out.println(progress * 200000 + " input lines have been read !");
			}
			br.close();
			fr.close();
			
			//if sb is not empty
			if (sb.toString().trim().length() > 0){
				st = splitRegex.split(sb.toString());
				sb = new StringBuilder();
				
				if (st.length % nbr != 0){
					throw new Exception("There is a line that doesn't contain exactly 3 numbers");
				}
				//each line contains 3 numbers
				for(int i = 0; i < st.length / nbr; i++){
					x = Integer.parseInt(st[i * nbr + 0]);
					
					//position2
					y = Integer.parseInt(st[i * nbr + 1]);
					
					//keeping absolute positions, so that later they can be recovered from indices
					
					posSet.add(x);
					
					posSet.add(y);
					 
				}
				
			}

			lstPos.addAll(posSet);
			
			//sort the lst of position ascendingly
			Collections.sort(lstPos);			
			
			//map postion into absolute index
			for(int i = 0; i < lstPos.size(); i++){
				hmPos.put(lstPos.get(i), i);
			}
			
			//initialize the matrix contact
			a = new double[lstPos.size()][lstPos.size()];			
			
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			count = 0;
			progress = 0;
			while((ln = br.readLine()) != null){
				
				if (ln.startsWith("#") || ln.trim().length() == 0){
					continue;
				}
				
				//read every a thoudsand lines and split at once
				sb.append(ln).append(" ");
				count++;
				
				if (count == 200000){
					count = 0;
					st = splitRegex.split(sb.toString());
					sb = new StringBuilder();
					
					if (st.length % nbr != 0){
						throw new Exception("There is a line that doesn't contain exactly 3 numbers");
					}
					//each line contains 3 numbers
					for(int i = 0; i < st.length / nbr; i++){
						x = Integer.parseInt(st[i * nbr + 0]);
						
						//position2
						y = Integer.parseInt(st[i * nbr + 1]);
						
						//interaction frequency
						f = Double.parseDouble(st[i * nbr + 2]);
						
						//keeping absolute positions, so that later they can be recovered from indices
						
						id1 = hmPos.get(x);
						id2 = hmPos.get(y);
												
						//if the frequency is less than thresholds, ignore the contact so that it is considered as non-contact
						

						a[id1][id2] = f;
						a[id2][id1] = f;
						
						 
					}
					
					progress++;
					//System.out.println(progress * 200000 + " input lines have been read - second time!");
				}	

			}

			//if sb is not empty
			if (sb.toString().trim().length() > 0){
				st = splitRegex.split(sb.toString());
				sb = new StringBuilder();
				
				if (st.length % nbr != 0){
					throw new Exception("There is a line that doesn't contain exactly 3 numbers");
				}
				//each line contains 3 numbers
				for(int i = 0; i < st.length / nbr; i++){
					x = Integer.parseInt(st[i * nbr + 0]);
					
					//position2
					y = Integer.parseInt(st[i * nbr + 1]);
					
					//interaction frequency
					f = Double.parseDouble(st[i * nbr + 2]);
					
					//keeping absolute positions, so that later they can be recovered from indices
					
					id1 = hmPos.get(x);
					id2 = hmPos.get(y);
					
					a[id1][id2] = f;
					a[id2][id1] = f;
					
					 
				}
				
			}

		
		}catch(FileNotFoundException ex){
			ex.printStackTrace();
			throw ex;
		}catch(IOException ex){
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
		
		return a;
	}

	
}
