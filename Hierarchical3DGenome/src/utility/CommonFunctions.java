package utility;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.List;

import valueObject.RegionVO;

public class CommonFunctions {

	
	public static void make_folder(String folder) throws Exception{
		File file = new File(folder);
		file.mkdirs();
	}
	
	public static boolean isExist(String fileName){
		File f = new File(fileName);
		return f.exists();
	}
	
	public static void delete_file(String file_name) throws Exception{
		File file = new File(file_name);
		if (file.exists() && file.isFile()) file.delete();
		else if(file.isDirectory()){
			for(File f : file.listFiles()){
				delete_file(f.getAbsolutePath());
			}
			file.delete();
		}
	}
	
	/**
	 * Given a folder, find the file contain the mapping between ID <-> genomic regions
	 * @param folder_name
	 * @return
	 */
	public static String find_mapping_file(String folder_name){
		File folder = new File(folder_name);
		if (!folder.exists()) return null;
		
		for(File f : folder.listFiles()){
			if (f.getName().contains("coordinate_mapping")){
				return f.getAbsolutePath();
			}
		}
		return null;
	}
	
	/**
	 * Given a folder, find the file contain the pdb file
	 * @param folder_name
	 * @return
	 */
	public static String find_pdb_file(String folder_name){
		File folder = new File(folder_name);
		if (!folder.exists()) return null;
		
		for(File f : folder.listFiles()){
			if (f.getName().endsWith(".pdb")){
				return f.getAbsolutePath();
			}
		}
		return null;
	}

	/**
	 * Given a folder, find the file contain the best convert factor
	 * @param folder_name
	 * @return
	 */
	public static String find_best_convert_factor_file(String folder_name){
		File folder = new File(folder_name);
		for(File f : folder.listFiles()){
			if (f.getName().contains("best_alpha")){
				return f.getAbsolutePath();
			}
		}
		return null;
	}
	

	/**
	 * Given a folder, find the file contain the best convert factor
	 * @param folder_name
	 * @return
	 */
	public static double read_best_convert_factor(String file_name) throws Exception{
		
		BufferedReader br = new BufferedReader(new FileReader(new File(file_name)));
		String ln, st[];
		while((ln = br.readLine()) != null){
			if (ln.length() > 0 && ln.startsWith("Best convert factor:")){
				st = ln.split("[:,]");
				br.close();
				return Double.parseDouble(st[1]);
			}
		}
		
		br.close();		
		return Double.NaN;
	}
	
	/**
	 * Get file name or last directory from a path
	 * @param path
	 * @return
	 */
	public static String get_file_name_from_path(String path){
		StringBuilder sb = new StringBuilder();
		for(int i = path.length() - 1; i >= 0; i--){
			if (path.charAt(i) == '/' || path.charAt(i) == '\\'){
				return sb.reverse().toString();
			}else{
				sb.append(path.charAt(i));
			}
		}
		
		return sb.toString();
	}
	
	/**
	 * Check if two regions overlap
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static int overlap(RegionVO p1, RegionVO p2){
		if (p1.getChr_id() != p2.getChr_id()) return 0;
		
		//overlap length
		int len = Math.min(p1.getEnd(), p2.getEnd()) - Math.max(p1.getStart(), p2.getStart()) + 1;
		
		return len;
		
//		if (len <= 0) return 0.0;
//		
//		//return the percentage of overlap 
//		return Math.min( len * 100.0 / (p1.getEnd() - p1.getStart()), len * 100.0 / (p2.getEnd() - p2.getStart()));
	}
	
	/**
	 * Look for a region in a sorted list
	 * @param lst
	 * @param r : object with start == end == genomic sequence id
	 * @return
	 */
	public static RegionVO lookup_region(List<RegionVO> lst, RegionVO r){
		
		//binary search at two levels
		int mid, i = 0, j = lst.size() - 1;
		RegionVO reg = null, mid_reg;
		while(i <= j){
			
			mid = i/2 + j/2 + (i % 2 + j % 2) / 2;
			mid_reg = lst.get(mid);
			if (overlap(mid_reg, r) > 0) return mid_reg;
			
			if (mid_reg.getChr_id() > r.getChr_id() || (mid_reg.getChr_id() == r.getChr_id() && mid_reg.getStart() > r.getEnd())){
				j = mid - 1;
			}else{
				i = mid + 1;				
			}			
		}
		
		return reg;
	}
	
}
