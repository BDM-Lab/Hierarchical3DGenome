package utility;

public class RegionVO implements Comparable<RegionVO> {
	
	int id; //id in relative with other regions, 
	int start; //starting id in its genomic sequence, inclusive
	int end; //ending id in its genomic sequence, inclusive
	int chr_id; //chromosome id
	boolean isDomain = false; //if the region is a TAD

	public RegionVO(int chr_id, int start, int end){		
		this.start = start;
		this.end = end;
		this.chr_id = chr_id;
	}	
	
	public RegionVO(int id, int chr_id, int start, int end){
		this.id = id;
		this.start = start;
		this.end = end;
		this.chr_id = chr_id;
	}
	
	public int compareTo(RegionVO vo){
		if (vo.chr_id != this.chr_id) return Integer.compare(this.chr_id, vo.chr_id);
		return Integer.compare(this.start, vo.start);
	}

	
	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public int getChr_id() {
		return chr_id;
	}

	public void setChr_id(int chr_id) {
		this.chr_id = chr_id;
	}

	public boolean isDomain() {
		return isDomain;
	}

	public void setDomain(boolean isD) {
		this.isDomain = isD;
	}		
}

