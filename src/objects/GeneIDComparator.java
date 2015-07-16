package objects;

import java.util.Comparator;

public class GeneIDComparator implements Comparator<GenomicCoordinate> {
	public int compare(GenomicCoordinate o1, GenomicCoordinate o2) {
		return o1.getAttribute("gene_id").compareTo(o2.getAttribute("gene_id"));
	}
}
