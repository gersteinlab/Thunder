package utils;

public class Object_utils {

	
	/**
	 * Inititalise array of Integer objects
	 * @param length size of the array
	 * @param initValue initial value placed in each cell of the array
	 * @return
	 */
	public static Integer[] initArray_Integer(int length, Integer initValue){
		Integer[] array = new Integer[length];
		for(int i=0;i<length;i++)
			array[i] = initValue;
		return array;
	}
	public static Integer[] initArray_Integer(int length){
		return initArray_Integer(length, 0);
	}
	
	/**
	 * Inititalise array of Double objects
	 * @param length size of the array
	 * @param initValue initial value placed in each cell of the array
	 * @return
	 */
	public static Double[] initArray_Double(int length, Double initValue){
		Double[] array = new Double[length];
		for(int i=0;i<length;i++)
			array[i] = initValue;
		return array;
	}
	public static Double[] initArray_Double(int length){
		return initArray_Double(length, 0.0);
	}
	
	
	/**
	 * Inititalise array of doubles
	 * @param length size of the array
	 * @param initValue initial value placed in each cell of the array
	 * @return
	 */
	public static double[] initArray_double(int length, double initValue){
		double[] array = new double[length];
		for(int i=0;i<length;i++)
			array[i] = initValue;
		return array;
	}
	public static double[] initArray_double(int length){
		return initArray_double(length, 0.0);
	}
	
	
	/**
	 * Inititalise array of ints
	 * @param length size of the array
	 * @param initValue initial value placed in each cell of the array
	 * @return
	 */
	public static int[] initArray_int(int length, int initValue){
		int[] array = new int[length];
		for(int i=0;i<length;i++)
			array[i] = initValue;
		return array;
	}
	public static int[] initArray_int(int length){
		return initArray_int(length, 0);
	}
	
	
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
