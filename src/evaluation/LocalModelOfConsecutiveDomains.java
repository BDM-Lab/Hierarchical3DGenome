package evaluation;

import java.util.ArrayList;
import java.util.List;

import full_modeling.HierarchicalModeling;
import noisy_mds.StructureGeneratorLorentz_HierarchicalModeling;
import topological_domain.IdentifyDomains;
import utility.CommonFunctions;
import valueObject.RegionVO;

/**
 * Generate models of regions consisting of 2 consecutive domains
 * @author Tuan
 *
 */
public class LocalModelOfConsecutiveDomains {

	public static void main(String[] args) throws Exception {
		
		String inputData = "/home/tuantrieu/lorentz/input/chr1_5kb_gm12878_list.txt";
		String domainFile = "/home/tuantrieu/lorentz/input/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_whole.txt";
		String domainInputFolder = "/home/tuantrieu/lorentz/output/chr1/large_domains_input";
		String domainOutputFolder = "/home/tuantrieu/lorentz/output/chr1/large_domains_output";

//		String inputData = "E:/GM12878/1kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_1kb_gm12878_list.txt";
//		String domainFile = "E:/GM12878/loop_arrow_domain/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_whole.txt";
//		String domainInputFolder = "output/chr1/large_domains_input";
//		String domainOutputFolder = "output/chr1/large_domains_output";

		int chr = 1;
		int res = 5000;
		
		if (args.length > 0){
			inputData = args[0];
			
			chr = Integer.parseInt(args[1]);
			res = Integer.parseInt(args[2]);
			
			domainFile = args[3];
			domainInputFolder = args[4];
			domainOutputFolder = args[5];
		}
		
				
		CommonFunctions.delete_file(domainInputFolder);
		CommonFunctions.make_folder(domainInputFolder);		
		
		IdentifyDomains identifyDomains = new IdentifyDomains(domainFile, chr, res);
		
		List<RegionVO> regions = identifyDomains.getDomains();
		regions = identifyDomains.merge2AdjacentDomains(regions);
		
		//System.out.println((regions.get(619).getEnd() - regions.get(619).getStart())/1000);
				
		HierarchicalModeling.extract_domain_contact(inputData, regions, domainInputFolder, chr, 0);
		
		double convert_factor = 0.9;// 0.9 for 5kb resolution, 1.2 for 1kb
		
		//build 3D models for domains
		String[] params = new String[2];
		params[0] = domainInputFolder;
		params[1] = domainOutputFolder;
		
		CommonFunctions.delete_file(domainOutputFolder);
		CommonFunctions.make_folder(domainOutputFolder);		
		StructureGeneratorLorentz_HierarchicalModeling.run_for_folder(params, convert_factor,1);
		
		CommonFunctions.delete_file(domainInputFolder);

	}

}
