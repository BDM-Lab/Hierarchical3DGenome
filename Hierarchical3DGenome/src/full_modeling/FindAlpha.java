package full_modeling;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.inference.TTest;

import topological_domain.IdentifyDomains;
import utility.RegionVO;

/**
 * Read converter factor to try to find a consensus value
 * @author Tuan
 *
 */
public class FindAlpha {
	
	public static void main(String[] args) throws Exception{
		String input_folder = "C:/Users/Tuan/workspace/GenomeMDS/output/hierarchical_modeling/chr1/backup";
		String input_domain_file = "E:/GM12878/loop_arrow_domain/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_whole.txt";
		int resolution = 1000;
		
		File folders = new File(input_folder);
		HashMap<Integer, Double> factor_domain = new HashMap<Integer, Double>();
		HashMap<Integer, Double> factor_nondomain = new HashMap<Integer, Double>();
		
		IdentifyDomains domain_identifer = new IdentifyDomains(input_domain_file, resolution);				
		List<RegionVO> regions = domain_identifer.get_all_regions();
		//to indicate if a region is a domain or non-domain
		HashMap<Integer, Boolean> region_domain = new HashMap<Integer, Boolean>();
		
		for(RegionVO reg : regions){
			region_domain.put(reg.getId(), reg.isDomain());
		}
		
		SummaryStatistics domain = new SummaryStatistics();
		SummaryStatistics nondomain = new SummaryStatistics();
		
		DescriptiveStatistics des_domain = new DescriptiveStatistics();
		DescriptiveStatistics des_nondomain = new DescriptiveStatistics();
		
		for(File folder : folders.listFiles()){
			for(File f : folder.listFiles()){
				if (f.getName().equals("best_alpha_log.txt")){
					BufferedReader br = new BufferedReader(new FileReader(f));
					String ln;
					String[] st;
					while((ln = br.readLine()) != null){
						if (ln.length() == 0) continue;
						st = ln.split(",");
						int regid = Integer.parseInt(folder.getName().replace("region_", ""));
						double factor = Double.parseDouble(st[0].replace("Best convert factor: ", ""));
						
						if (Math.abs(-1.0 - factor) < 1e-6){
							System.out.println(folder.getName());
							continue;
						}
						
						des_domain.addValue(factor);
						
						if (region_domain.get(regid)){
							domain.addValue(factor);
							
							
							factor_domain.put(regid,factor);							
						}else{
							factor_nondomain.put(regid,factor);
							des_nondomain.addValue(factor);
							nondomain.addValue(factor);
						}
					}
					br.close();
					break;
				}
			}
		}
		
		TTest ttest = new TTest();
		double pvalue = ttest.tTest(domain,  nondomain);
		
		System.out.println(pvalue);
		
		System.out.println(domain.getSummary());
		System.out.println(des_domain.getPercentile(50.0));
		
		System.out.println(nondomain.getSummary());
		System.out.println(des_nondomain.getPercentile(50.0));
		
		System.out.println(domain.getMean() - nondomain.getMean());
		

	}

}
