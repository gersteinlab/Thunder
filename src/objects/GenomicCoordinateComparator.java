package objects;

import java.util.Comparator;

public class GenomicCoordinateComparator implements Comparator<GenomicCoordinate> {
	public int compare(GenomicCoordinate o1, GenomicCoordinate o2) {
		return o1.getStart().compareTo(o2.getStart());
	}
}
