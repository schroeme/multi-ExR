// ************* FIJI SCRIPT FOR MEASURING AVERAGE INTENSITY OF A MIP ************* \\
// Last modified by MES on 11/27/23
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

parent = "E:/Margaret/mExR/2022.05_validation/zcropped/"; //source directory
setBatchMode(true); //true will increase speed
run("Close All");
image_folder = parent;

list = getFileList(image_folder);									 				
extension = "tif";	
contains = "ch04";
index=-1;												
    	
for (i=0; i<list.length; i++) {								 //go through all files in the current data folder
	if (endsWith(toLowerCase(list[i]), "." + extension)){			//if its extension matches up with the extensions of interest
		if (indexOf(list[i],contains) >= 0){
			
			open(image_folder+list[i]);	
			run("Z Project...", "projection=[Max Intensity]");				
			//print(title);
			run("Measure");
			close();
		} // end if loop selecting only files with "contains" string
	} // end "if" loop selecting only files with the right extension
} //end "for" loop going through all files in current folder 

