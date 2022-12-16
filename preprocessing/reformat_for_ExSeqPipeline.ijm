// ************* FIJI SCRIPT FOR REFORMATTING PREPROCESSED .TIF STACKS FOR EXSEQPROCESSING REGISTRATION ************* \\
// Last modified by MES on 12/15/22
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

setBatchMode(true); // for speed
run("Close All");

//source directory
parentfolder = "/run/user/1000/gvfs/smb-share:server=exr.mit.edu,share=exr/Margaret/mExR/2022.10_5xFAD/preprocessed/"; //CHANGE HERE
//target directory
targetfolder = "/ssd/mExR/2022.10_5xFAD/0_raw/"; //CHANGE HERE. Please note the folder structure requred by ExSeqProcessing

image_folder = parentfolder;
list = getFileList(image_folder);
print(image_folder);//find all files in current data folder	
extension = "tif";							 //set file types to look at 'tif' stacks
index=-1; 									//initialize index that will mark files with 
	    	
print(list.length);
for (i=0; i<list.length; i++) {								 //go through all files in the current data folder
	if (endsWith(toLowerCase(list[i]), "." + extension)){			//if its extension matches up with the extensions of interest
		//open tiff stack
		open(image_folder + list[i]);
		index=i;
		name =  list[index];												//record the full file name (including file extension; eg:'stack1.tif')
        dotIndex = lastIndexOf(name, ".");										//record position of dot that separates extension from filename	
        title = substring(name, 0, dotIndex); 								//record the abbreviated file name (excluding file extension; eg:'stack1')
        print(title);
        splits = split(title, "_");
        basename=splits[0]; //may have to change this depending on dataset
        
        getDimensions(w, h, channels, slices, frames);
        
    	roundname = "_" + splits[1];
    	
        run("Split Channels");

        if (channels==4) {
	        saveAs("Tiff", targetfolder + basename + roundname + "_ch04.tif");
	        close(); 
	        //saveAs("Tiff", targetfolder + basename + roundname + "_ch03.tif"); 
	        //close();
	        saveAs("Tiff", targetfolder + basename + roundname + "_ch03.tif"); 
	        close();
	        saveAs("Tiff", targetfolder + basename + roundname + "_ch02.tif"); 
	        close();
	        saveAs("Tiff", targetfolder + basename + roundname + "_ch01.tif"); 
	        close();
        } else { //CHANGE HERE if you have an uneven number of channels between rounds
        	//Currently, this is set up to duplicate one of the channels
        	
        	//saveAs("Tiff", targetfolder + basename + roundname + "_ch04.tif");
	        saveAs("Tiff", targetfolder + basename + roundname + "_ch03.tif"); 
	        close();
	        saveAs("Tiff", targetfolder + basename + roundname + "_ch02.tif"); 
	        saveAs("Tiff", targetfolder + basename + roundname + "_ch01.tif"); 
	        close();
	        saveAs("Tiff", targetfolder + basename + roundname + "_ch04.tif"); 
	        close();
        }
        
	    run("Close All");
	    run("Collect Garbage");
		
	} // end "if" loop selecting only files with the right extension
} //end "for" loop going through all files in current folder 


run("Clear Results");
