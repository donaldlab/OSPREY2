
public class Util {
	// Returns the min of an array
	public static double min(double array []){
		int minIndex = 0;
		for (int i = 1; i < array.length; i++){
			if(array[minIndex] > array[i]){
				minIndex = i;
			}			
		}
		return array[minIndex];
	}
}
