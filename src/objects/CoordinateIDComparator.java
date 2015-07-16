package objects;

import java.util.Comparator;

public class CoordinateIDComparator implements Comparator<GenomicCoordinate> {
	public int compare(GenomicCoordinate o1, GenomicCoordinate o2) {
		return o1.getCoordinateID().compareTo(o2.getCoordinateID());
	}
}
