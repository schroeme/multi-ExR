// ************* FIJI SCRIPT FOR SHORTENING MULTI-EXR Z-STACKS AFTER REGISTRATION ************* \\
// Last modified by MES on 12/15/22
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

parent = "E:/Margaret/mExR/2022.10_5xFAD/registered_v2/"; //source directory
targetdir = "E:/Margaret/mExR/2022.10_5xFAD/cropped_z/"; //target directory
setBatchMode(true); //true will increase speed
run("Close All");
image_folder = parent;

list = getFileList(image_folder);									 				
extension = "tif";		
index=-1;												
    	
for (i=0; i<list.length; i++) {								 //go through all files in the current data folder
	if (endsWith(toLowerCase(list[i]), "." + extension)){			//if its extension matches up with the extensions of interest
			open(image_folder+list[i]);
			index=i;
			name =  list[index];												
	        dotIndex = indexOf(name, ".");										
	        title = substring(name, 0, dotIndex); 								
	        
			range="10-70"; //CHANGE HERE for desired z-range. 
			print(title);
			run("Duplicate...", "duplicate range=" + range + " use");
			
			saveAs("Tiff", targetdir + title + ".tif"); // save the channel we want to extract from separately
			close();
	
		    run("Close All");
		    run("Collect Garbage");
		//}
	} // end "if" loop selecting only files with the right extension
} //end "for" loop going through all files in current folder 

		
run("Clear Results");
