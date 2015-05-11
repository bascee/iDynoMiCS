package utils;

public final class Helpers {


	public static float[] doublesToFloats(double[] dArray) {
		float[] fArray = new float[dArray.length];
		for (int i=0; i<dArray.length; i++) {
			fArray[i] = (float) dArray[i];
		}
		return fArray;
	}
	
	public static float[] fillArray(float filling, int length) {
		float[] fArray = new float[length];
		for (int i=0; i<length; i++) {
			fArray[i] = filling;
		}
		return fArray;
	}
}
