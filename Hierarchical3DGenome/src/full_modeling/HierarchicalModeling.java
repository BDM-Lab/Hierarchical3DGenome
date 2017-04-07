package full_modeling;

import static utility.CommonFunctions.lookup_region;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import noisy_mds.StructureGeneratorLorentz_HierarchicalModeling;
import topological_domain.IdentifyDomains;
import utility.CommonFunctions;
import utility.DistanceEstimate;
import utility.ExtractWithinAndConsecutiveDomainContact;
import utility.Helper;
import utility.RegionVO;
import valueObject.Constants;
import valueObject.Constraint;
import valueObject.InputParameters;

public class HierarchicalModeling {
	
	//private static String input_folder = "E:/GM12878/1kb_resolution_intrachromosomal/chr1/MAPQGE30/";
	//private static String input_data_file = input_folder + "/chr1_1kb.RAWobserved";	
	//private static String input_data_domain_file = input_folder + "/chr1_1kb_gm12878_list_domain.txt";
	
	private static Helper helper = Helper.getHelperInstance();
	private static int resolution = 1000;	
	private static int chr_id = 2;
	
	public static void main(String[] args) throws Exception{
		
		String input_data_file, input_observed_data_file, input_domain_file, output_folder;
		if (args.length == 0){
			input_data_file = "E:/GM12878/1kb_resolution_intrachromosomal/chr2/MAPQGE30/chr2_1kb_gm12878_list.txt";
			input_observed_data_file = "E:/GM12878/1kb_resolution_intrachromosomal/chr2/MAPQGE30/chr2_1kb.RAWobserved";
			input_domain_file = "E:/GM12878/loop_arrow_domain/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_whole.txt";
			output_folder = "output/test/";
		}else{
			chr_id = Integer.parseInt(args[0]);
			input_data_file = args[1];
			input_observed_data_file = args[2];
			input_domain_file = args[3];
			output_folder = args[4];
		}		
		
		File file = new File(output_folder);
		if (!file.exists()) file.mkdirs();
		
		
		String output_folder_low_res_global = output_folder + "/low_res_global/";
		
		String file_name = CommonFunctions.get_file_name_from_path(input_data_file);
		
		file_name = file_name.contains(".") ? file_name.split("\\.")[0] : file_name;
		
		//contact file whose points are domains, 
		String unnorm_domain_contact_file = output_folder + "/" + file_name + "_domain_contact_unnorm.txt";
		String domain_contact_file = output_folder + "/" + file_name + "_domain_contact.txt";
		
		//new contact data file with contacts between far away domains removed
		String filtered_contact_data_file = output_folder + "/" + file_name + "_filtered.txt";
		
		//String filtered_observed_contact_data_file = output_folder + "/observed_" + file_name + "_filtered.txt";
		
		//folder contains within-domain contacts, each file is for one domain
		String contact_domain_file_input_folder = output_folder + "/domains_input/";
		
		//folder contains 3D models of domains
		String contact_domain_file_output_folder = output_folder + "/domains_output/";
		
		//folder to contain global 3D model
		
		compute_contact_map_of_domains(input_observed_data_file, input_domain_file, unnorm_domain_contact_file);
		//normalize the contact domain file
		normalization.Normalizer.normalize(unnorm_domain_contact_file, domain_contact_file, 1);//normalize with ICE
				
		CommonFunctions.delete_file(contact_domain_file_input_folder);
		CommonFunctions.make_folder(contact_domain_file_input_folder);
		extract_domain_contact(input_data_file, input_domain_file, contact_domain_file_input_folder);
		
		
		//generate global model at low resolution (each point is a domain)
		InputParameters input_parameter = new InputParameters();			
		input_parameter.setNum(1);				
		input_parameter.setOutput_folder(output_folder_low_res_global);		
		input_parameter.setInput_file(domain_contact_file);		
		input_parameter.setFile_prefix(file_name);
		
		input_parameter.setVerbose(false);
		input_parameter.setLearning_rate(1.0);
		input_parameter.setMax_iteration(20000);		
		input_parameter.setKeepOriginalScale(true);
		
		input_parameter.setNumber_threads(10); 
		
		//input_parameter.setConvert_factor(0.6);//for chr.1		
		CommonFunctions.delete_file(output_folder_low_res_global);
		StructureGeneratorLorentz_HierarchicalModeling.run_for_one_input(input_parameter);
		//
		
		//get model and mapping files
		String global_model_file = CommonFunctions.find_pdb_file(output_folder_low_res_global);
		String global_mapping_file = CommonFunctions.find_mapping_file(output_folder_low_res_global);
		//
		
		
		//filter contacts from 2 far away domains
		ExtractWithinAndConsecutiveDomainContact.determine_proximal_domain_and_filter(input_data_file, input_domain_file, 
								global_model_file, global_mapping_file, filtered_contact_data_file, chr_id);
		
		//now, I can delete the original input_data_file to have more space		
		

//		System.out.println("Ratio: " + ratio);
		
		double convert_factor = 1.0, ratio = 1.0;
		
		
		//PrintWriter pw = new PrintWriter(output_folder + "/convertFactorSearch.txt");
		
		/*
		int trial = 20;			
		ArrayList<Double> cf_lst = new ArrayList<Double>();		
		for(int i = 0; i < trial; i++){
			
			int[] range = determine_chunk_size(input_domain_file,2);
			
			System.out.println("Chunk range:" + range[0] + "\t" + range[1]);
			//pw.println("Chunk range:" + range[0] + "\t" + range[1]);
			
//			range[0] = 116115001;         
//			range[1] = 117899999; 
			//Extract a chunk of data to find the best convert_factor alpha
			String contact_data_chunk = output_folder + "/" + file_name + "_chunk.txt";
			
			extract_contact_data(filtered_contact_data_file, contact_data_chunk, range[0], range[1]);						
			
			//now, build 3D model with this chunk of data and look for best convert_factor alpha
			String chunk_data_output_folder = output_folder + "/chunk_data/";
			input_parameter = new InputParameters();
			input_parameter.setNum(1);				
			input_parameter.setOutput_folder(chunk_data_output_folder);		
			input_parameter.setInput_file(contact_data_chunk);		
			input_parameter.setFile_prefix(CommonFunctions.get_file_name_from_path(contact_data_chunk).split("\\.")[0]);
			
			input_parameter.setVerbose(false);
			input_parameter.setLearning_rate(0.1);
			input_parameter.setMax_iteration(20000);		
			input_parameter.setKeepOriginalScale(true);
			input_parameter.setAddInequalityConstraint(false);
	
			//input_parameter.setConvert_factor(1.2);			
			//input_parameter.setNumber_threads(1);//only in local
			
			CommonFunctions.delete_file(chunk_data_output_folder);			
			double cor = StructureGeneratorLorentz_HierarchicalModeling.run_for_one_input(input_parameter);
						
			convert_factor = CommonFunctions.read_best_convert_factor(CommonFunctions.find_best_convert_factor_file(chunk_data_output_folder));
			cf_lst.add(convert_factor);
			
			pw.println("convert_factor:" + convert_factor + "\tcor:" + cor);
			System.out.println("convert_factor:" + convert_factor + "\tcor:" + cor);
		}
		
		Collections.sort(cf_lst);
		pw.println();
		for(int i = 0; i < cf_lst.size(); i++){
			pw.printf("%.2f, ", cf_lst.get(i));
		}
		
		pw.println();
		
		convert_factor = cf_lst.get((cf_lst.size() - 1)/ 2);
		
		convert_factor = 1.2;
		
		pw.println();
		pw.println("Convert factor:" + convert_factor);
		System.out.println("Convert factor:" + convert_factor);
		pw.flush();
		
		*/
		
		convert_factor = 1.2;
		
		//build 3D models for domains
		String[] params = new String[2];
		params[0] = contact_domain_file_input_folder;
		params[1] = contact_domain_file_output_folder;
		
		CommonFunctions.delete_file(contact_domain_file_output_folder);
		CommonFunctions.make_folder(contact_domain_file_output_folder);		
		StructureGeneratorLorentz_HierarchicalModeling.run_for_folder(params, convert_factor);
		
		
		
		//calculate the ratio
		ArrayList<Integer> lstPos = new ArrayList<Integer>();
		List<Constraint> lstCons = helper.readContactList(filtered_contact_data_file, lstPos,0.0);
		
		IdentifyDomains identifyDomain = new IdentifyDomains(input_domain_file, chr_id, resolution);
		identifyDomain.assignDomainID(lstCons, identifyDomain.getDomains());
		
		List<Constraint> lstConsDomain = identifyDomain.extractContactBetweenDomain(lstCons, 10,11);
		
		DistanceEstimate.distanceBewteen2Domains(lstConsDomain, convert_factor, 10, 11, contact_domain_file_output_folder);
		
		//pw.println();
		ratio = DistanceEstimate.distanceEstimateBetweenDomains(lstCons, convert_factor, global_model_file, 
				input_domain_file, contact_domain_file_output_folder, chr_id, resolution);
		
		//pw.println("Ratio:" + ratio);	
		//pw.close();
		
		//ratio = 2.569;
		
		
		//	
		
		init_structure(filtered_contact_data_file, input_domain_file, global_model_file, global_mapping_file, contact_domain_file_output_folder, convert_factor, ratio, output_folder);
		
		CommonFunctions.delete_file(filtered_contact_data_file);
		CommonFunctions.delete_file(input_data_file);
		CommonFunctions.delete_file(unnorm_domain_contact_file);
		
	}
	
	
	
	
	public static void init_structure(String input_file, String input_domain_file, String global_model_file, String global_coordinate_mapping_file, 
			String loci_model_folder, double convert_factor, double global_model_scale, String output_folder) throws Exception{
		
		
		IdentifyDomains domain_identifer = new IdentifyDomains(input_domain_file, chr_id, resolution);		
		List<RegionVO> regions = domain_identifer.getDomains();		
		
		Helper helper = Helper.getHelperInstance();
		
		ArrayList<Integer> lstPos = new ArrayList<Integer>();
		List<Constraint> lstCons = helper.readContactList(input_file, lstPos,0.0);
		HashMap<Integer, Integer> map_id_cor = new HashMap<Integer, Integer>();
		for(int i = 0; i < lstPos.size(); i++){
			map_id_cor.put(i, lstPos.get(i)); // mapping id to coordinate
		}
		lstPos = null; // marked memory is unused
		
		double scale_factor_for_local_model = 1.0;// is the average of converted distances
		int number_points = map_id_cor.size();//compute_string_length(input_file);//219899 - full file
		//int number_points = 20000;
		
		//calculate scale_factor = 0.8359535263264956  for chr.1	, 1.0377017447829227 for first 20000
		scale_factor_for_local_model = compute_scale_factor(lstCons, convert_factor);
		
		System.out.println("scale_factor_for_local_model:" + scale_factor_for_local_model);
		
		double[][] global_model = helper.loadPDBStructure(global_model_file);		
		
		System.out.println("Global model scale:" + global_model_scale);
		
		for(int i = 0; i < global_model.length; i++){
			global_model[i][0] *= global_model_scale;
			global_model[i][1] *= global_model_scale;
			global_model[i][2] *= global_model_scale;
		}
		//
		
		//chr1_1kb_gm12878_list_domain_norm_coordinate_mapping.txt
		BufferedReader br = new BufferedReader(new FileReader(new File(global_coordinate_mapping_file)));
		String ln;
		String[] st;
		HashMap<Integer, Integer> map_domain_id_cor = new HashMap<Integer, Integer>();
		while((ln = br.readLine()) != null){
			st = ln.split("\\s+");
			map_domain_id_cor.put(Integer.parseInt(st[0]), Integer.parseInt(st[1]));
		}
		br.close();
		//
		
		double[] str = new double[number_points * 3];
		
		int current = 0;
		
		int local_skip = 0;
		String folder_name;
		double loci_model[][], prevX=0, prevY=0, prevZ=0;
		String mapping_file;
		HashMap<Integer,Integer> map_id_cor_tmp; //contains mapping of indices to coordinates
		for(int t = 0; t < regions.size(); t++){
			folder_name = loci_model_folder + "/region_" + regions.get(t).getId();
			File folder = new File(folder_name);
			
			loci_model = null;
			
			if (!folder.exists()) continue;
			
			mapping_file = CommonFunctions.find_mapping_file(folder_name);
			
			ArrayList<String> fileNames = new ArrayList<String>();
			for(File f : folder.listFiles()){				
				if (f.getName().endsWith(".pdb")){
					fileNames.add(f.getAbsolutePath());
				}
			}
			if (fileNames.size() == 0){
				System.out.println("No model in: " + folder_name);
				continue;
			}
			
			String modelFile = fileNames.get( (int) (Math.random() * fileNames.size()));
				
				
					
					loci_model = helper.loadPDBStructure(modelFile);					
					
					map_id_cor_tmp = helper.loadCoordinateMapping(mapping_file);
					
					//System.out.println(f.getName() + "\t" + map_id_cor_tmp.get(0) + "\t" + map_id_cor_tmp.get(Collections.max(map_id_cor_tmp.keySet())));
					
					//scale local models to the new scale of the whole model
					loci_model = helper.zoomStructure(loci_model, 1.0/scale_factor_for_local_model);
					
					double cx=0, cy=0, cz=0;
					for(int i = 0; i < loci_model.length; i++){
						cx += loci_model[i][0];
						cy += loci_model[i][1];
						cz += loci_model[i][2];
					}
					
					cx /= loci_model.length;
					cy /= loci_model.length;
					cz /= loci_model.length;
					
					int global_id = map_domain_id_cor.get(regions.get(t).getId());
					
					double tx = global_model[global_id][0] - cx;
					double ty = global_model[global_id][1] - cy;
					double tz = global_model[global_id][2] - cz;
					
					for(int i = 0; i < loci_model.length; i++){
						
						//if (current >= number_points) break;
						if (map_id_cor.get(current).intValue() == map_id_cor_tmp.get(i).intValue()){
							str[current * 3 + 0] = loci_model[i][0] + tx;
							str[current * 3 + 1] = loci_model[i][1] + ty;
							str[current * 3 + 2] = loci_model[i][2] + tz;
							prevX = str[current * 3 + 0];
							prevY = str[current * 3 + 1];
							prevZ = str[current * 3 + 2];
							current++;
							
						}else if(map_id_cor.get(current) < map_id_cor_tmp.get(i)){ // local model skip a point, increase global model index by 1
							
							while (map_id_cor.get(current) < map_id_cor_tmp.get(i)) {
								str[current * 3 + 0] = prevX + (Math.random() - 0.5);
								str[current * 3 + 1] = prevY + (Math.random() - 0.5);
								str[current * 3 + 2] = prevZ + (Math.random() - 0.5);
								
								prevX = str[current * 3 + 0];
								prevY = str[current * 3 + 1];
								prevZ = str[current * 3 + 2];
								
								
								//System.out.println("local model skips a point:" + current);
								current++;
								local_skip++;
							}
							
							if (map_id_cor.get(current).intValue() == map_id_cor_tmp.get(i).intValue()){
								str[current * 3 + 0] = loci_model[i][0] + tx;
								str[current * 3 + 1] = loci_model[i][1] + ty;
								str[current * 3 + 2] = loci_model[i][2] + tz;
								prevX = str[current * 3 + 0];
								prevY = str[current * 3 + 1];
								prevZ = str[current * 3 + 2];
								
								current++;
							}
							
							
						}else{
							System.err.println("global model skip a point but local model doesn't, error!:" + map_id_cor.get(current) + "\t" + map_id_cor_tmp.get(i));
							System.exit(1);
						}
						
					}					
					
				//	break; // to pick only one file
				//}
				//if (current >= number_points) break;
			//}					
			
		}
		
		System.out.println("Total of skipped points in local models:" + local_skip);
		
		if (current < number_points){
			System.err.println("Initial model is smaller than expected:" + current + "\t" + number_points);
			//System.exit(1);
		}
		
		for(int i = current; i < number_points; i++){			
			str[i * 3 + 0] = prevX + (Math.random() - 0.5);
			str[i * 3 + 1] = prevY + (Math.random() - 0.5);
			str[i * 3 + 2] = prevZ + (Math.random() - 0.5);
			
			prevX = str[i * 3 + 0];
			prevY = str[i * 3 + 1];
			prevZ = str[i * 3 + 2];			
		}
		
		//done with using contact for now
		//lstCons = null;
		
		//check if any distance is zero, make it non-zero for optimization 
		for(int i = 0; i < number_points; i++){
			for(int j = i + 1; j < number_points; j++){
				double d = helper.calEuclidianDist(str[i * 3], str[i * 3 + 1], str[i * 3 + 2], str[j * 3], str[j * 3 + 1], str[j * 3 + 2]);
				if (Math.abs(d) < 0.0001){
					str[j * 3] += Math.random() - 0.5;
					str[j * 3 + 1] += Math.random() - 0.5;
					str[j * 3 + 2] += Math.random() - 0.5;
				}
			}
		}
		//
		
		
		
		String init_file = output_folder + "/" + chr_id + "_init.pdb";
		
		helper.writeStructure(init_file, str, null, "test", false);		
		
		
		String tmpFolder = output_folder + "/tmp_model/";
		File file = new File(tmpFolder);
		file.mkdirs();
		
		
		
		InputParameters input_parameter = new InputParameters();
		
		input_parameter.setNum(10);	
		
		input_parameter.setOutput_folder(output_folder);
		
		input_parameter.setInput_file(input_file);
		
		//input_parameter.setContact_thres(1.0);
		
		input_parameter.setFile_prefix("chr" + chr_id);
		input_parameter.setConvert_factor(convert_factor);
		input_parameter.setVerbose(true);
		input_parameter.setLearning_rate(0.1);
		input_parameter.setMax_iteration(50000);
		
		input_parameter.setInitial_str_file(init_file);
		
		input_parameter.setAddInequalityConstraint(true);
		
		//input_parameter.setNumber_threads(1);
		
		//input_parameter.setTmpFolder(tmpFolder); //set to write out models during optimization
		
		StructureGeneratorLorentz_HierarchicalModeling generator = new StructureGeneratorLorentz_HierarchicalModeling(input_parameter);		
		generator.generateStructure();
		
		
	}
	
	/**
	 * 
	 * @param input_domain_file
	 * @return
	 */
	public static int[] determine_chunk_size(String input_domain_file, int ... gap) throws Exception{
		IdentifyDomains domain_identifer = new IdentifyDomains(input_domain_file, chr_id, resolution);		
		List<RegionVO> regions = domain_identifer.get_all_regions();
		int[] range = new int[2];
		
		int g = 0;
		if (gap.length > 0) g = gap[0];
		else g = 5;
		
		while(true){
			int id = (int) (Math.random() * (regions.size() - 5));
			//chunk span 2 adjacent regions
			range[0] = (regions.get(id).getStart() + regions.get(id).getEnd())/2;			
			range[1] = (regions.get(id + 1).getStart() + regions.get(id + 1).getEnd())/2;
			
			if ((range[1] - range[0])/resolution < 1000) break;
		
		}
		
		return range;
		
	}	
	
	/**
	 * Compute average distance from all contacts
	 * @param input_file
	 * @param convert_factor
	 * @return
	 * @throws Exception
	 */
	public static double compute_scale_factor(String input_file, double convert_factor) throws Exception{
		double freq,scale_factor = 0.0;
		int count = 0;
		
		BufferedReader br = new BufferedReader(new FileReader(new File(input_file)));
		String ln;
		String[] st;		
		while((ln = br.readLine()) != null){
			st = ln.split("\\s+");
			freq = Double.parseDouble(st[2]);
			
			scale_factor += (1.0 / Math.pow(freq,convert_factor));
			count++;
		}
		br.close();
		//		
		return scale_factor / count;
	}
	
	/**
	 * Compute average distance from all contacts
	 * @param lstCons
	 * @param convert_factor
	 * @return
	 * @throws Exception
	 */
	public static double compute_scale_factor(List<Constraint> lstCons, double convert_factor) throws Exception{
		double scale_factor = 0.0;
		
		for(Constraint con : lstCons){			
			scale_factor += (1.0 / Math.pow(con.getIF(),convert_factor));			
		}				
		return scale_factor / lstCons.size();
	}
	
	
	/**
	 * Compute the structure length to initialize the model
	 * @param input_file
	 * @param convert_factor
	 * @return
	 * @throws Exception
	 */
	public static int compute_string_length(String input_file) throws Exception{
		
		HashSet<Integer> set = new HashSet<Integer>();
		
		BufferedReader br = new BufferedReader(new FileReader(new File(input_file)));
		String ln;
		String[] st;	
		int x, y;
		while((ln = br.readLine()) != null){
			st = ln.split("\\s+");
			x = Integer.parseInt(st[0])/resolution;
			y = Integer.parseInt(st[1])/resolution;
			
			if (x == y) continue;
			
			set.add(x);
			set.add(y);
		}
		br.close();
		//		
		return set.size();
	}
	
	
	/**
	 * 1. Compute the contact matrix between domains
	 * 
	 * @param input_data_file
	 * @param input_domain_file
	 * @return
	 * @throws Exception
	 */
	public static void compute_contact_map_of_domains(String input_data_file, String input_domain_file, String domain_contact_file) throws Exception{
		
		IdentifyDomains domain_identifer = new IdentifyDomains(input_domain_file, chr_id, resolution);
		
		List<RegionVO> regions = domain_identifer.getDomains();
		Collections.sort(regions);

		//
		
		List<Constraint> constraint_list = new ArrayList<Constraint>();
		BufferedReader br = new BufferedReader(new FileReader(new File(input_data_file)));
		String ln;
		String[] st;
		RegionVO tmp = new RegionVO(chr_id,0,0), domain1, domain2;
		Constraint tmp_constraint = new Constraint(0,0,0);
		int x, y, i;
		double f;

		while((ln = br.readLine()) != null) {
			st = ln.split("\\s+");
			x = Integer.parseInt(st[0]);
			y = Integer.parseInt(st[1]);
			f = Double.parseDouble(st[2]);
			
			tmp.setChr_id(chr_id);
			
			//look for domain that contains first index
			tmp.setStart(x);
			tmp.setEnd(x);			
			domain1= lookup_region(regions, tmp);
			
			//look for domain that contains second index
			tmp.setStart(y);
			tmp.setEnd(y);
			domain2 = lookup_region(regions, tmp);
			
			tmp_constraint.setPos1(domain1.getId());
			tmp_constraint.setPos2(domain2.getId());
						
			i = Collections.binarySearch(constraint_list, tmp_constraint);
			
			
			if (i >= 0) constraint_list.get(i).setIF(constraint_list.get(i).getIF() + f);
			else {
				constraint_list.add(-i - 1, new Constraint(domain1.getId(),domain2.getId(), f));
			}
			
		}
		//pw.close();
		br.close();
		
		//output data
		PrintWriter pw = new PrintWriter(domain_contact_file);
		for(Constraint con : constraint_list){
			pw.println(con.getPos1() + "\t" + con.getPos2() + "\t" + con.getIF());
		}
		pw.close();		
	}

	/**
	 * 
	 * extract and print the intra-chromosomal contact matrix of each domain
	 * 
	 * @param input_data_file
	 * @param input_domain_file
	 * @return
	 * @throws Exception
	 */
	public static void extract_domain_contact(String input_data_file, String input_domain_file, String output_folder) throws Exception{
		
		IdentifyDomains domain_identifer = new IdentifyDomains(input_domain_file, chr_id, resolution);
		
		List<RegionVO> regions = domain_identifer.get_all_regions();
		Collections.sort(regions);
		
		//output files for contact domain of each region
		PrintWriter[] pws = new PrintWriter[regions.get(regions.size() - 1).getId() + 1];
		File tmp_folder = new File(output_folder);
		if (!tmp_folder.exists()) tmp_folder.mkdir();
		for(int i = 0; i < regions.size(); i++){
			//if (regions.get(i).isDomain()){
				pws[regions.get(i).getId()] = new PrintWriter(output_folder + "/" + "region_" + regions.get(i).getId() + ".txt");
			//}else{
				//pws[regions.get(i).getId()] = new PrintWriter(input_folder + "/tmp/" + "region_" + regions.get(i).getId() + ".txt");
			//}
		}
		//
		
		
		BufferedReader br = new BufferedReader(new FileReader(new File(input_data_file)));
		String ln;
		String[] st;
		RegionVO tmp = new RegionVO(1,0,0), domain1, domain2;
		
		int x, y, i;
		while((ln = br.readLine()) != null) {
			st = ln.split("\\s+");
			x = Integer.parseInt(st[0]);
			y = Integer.parseInt(st[1]);
			
			
			tmp.setChr_id(chr_id);
			
			//look for domain that contains first index
			tmp.setStart(x);
			tmp.setEnd(x);			
			domain1= lookup_region(regions, tmp);
			
			//look for domain that contains second index
			tmp.setStart(y);
			tmp.setEnd(y);
			domain2 = lookup_region(regions, tmp);
			
			//a contact is in a domain (not across domains)
			if (domain1.getId() == domain2.getId()){
				pws[domain1.getId()].println(ln);
			}			
			
		}
		//pw.close();
		br.close();
		
		for(i = 0; i < regions.size(); i++){
			pws[regions.get(i).getId()].close();;
		}		
		
	}
	
	
	
	public static void extract_contact_data(String input_file, String output_file, int start, int end) throws Exception{
		
		PrintWriter pw = new PrintWriter(output_file);
		BufferedReader br = new BufferedReader(new FileReader(new File(input_file)));
		String ln;
		String[] st;
		int x, y;
		long count = 0;
		while((ln = br.readLine()) != null){
			st = ln.split("\\s+");
			x = Integer.parseInt(st[0]);
			y = Integer.parseInt(st[1]);			
			
			if (x >= start && x <= end && y >= start && y <= end){
				pw.println(ln);
				count++;
				
				if (count % 100000 == 0) pw.flush();
			}
		}
		br.close();
		pw.close();		
		
	}
	

}
