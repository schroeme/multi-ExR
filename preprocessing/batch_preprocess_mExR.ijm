// ************* FIJI SCRIPT FOR PRE-PROCESSING MULTI-EXR IMAGES ************* \\
// Last modified by MES on 12/15/22
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

run("Close All");

parent = "E:/Jinyoung/2022.9.30 multi-ExR 5xFAD WT cortex/"; //parent folder where raw microscope images are stored (e.g. in .nd2 format)
parentout = "E:/Margaret/mExR/2022.10_5xFAD/"; //where to save the images
reg = "ctx"; //specify the brain region here (you can delete this from the final filename in the last section if you want)
condList = newArray( //condition folders - you will need to modify this based on your folder structure
	"5xFAD/",
	"WT/"
	);

roundList = newArray( //round folders - again, will need to be modified depending on how you saved your data
	"R1/before +SMI/",
	"R2/",
	"R3/",
	"R4/"
	"R5/"
	"R6/",
	"R7/"
	"R8/"
	"R9/",
	"R10/",
	"R11/"
	);

setBatchMode(true); //true for speed

for (i = 0; i < condList.length; i++){
    for (j = 0; j < roundList.length; j++){
    	
    	folder = parent+roundList[j]+condList[i];
    	imgList = getFileList(folder);
    	
    	condsplits = split(condList[i],"/");
    	cond = condsplits[0];
    	//regsplits = split(condsplits[1],"/"); // uncomment this if you have region information in the directory structure, which will allow you to forgo hard-coding

		//CHANGE HERE per the naming convention of your files. This assumes "ROI[X] 40x.tif"
    	roundind = indexOf(roundList[j],"R");
    	slashind = indexOf(roundList[j],"/");
    	roundno = substring(roundList[j], roundind+1, slashind);
    	
    	for (k = 0; k < imgList.length; k++){
	    	if (endsWith(imgList[k], ".nd2")){ //CHANGE HERE for your filetype
	    		if (indexOf(imgList[k],"40") >= 0){ //CHANGE HERE: provide a string that is ONLY in the images you're interested in
	    			//we have "40" because all of our high-mag images are named with "40x"
	    			// you can also skip this extra for-loop
	    			preprocessImage(imgList[k],folder,cond,reg,roundno);
				}
			}
    	}
	}
}


function preprocessImage(imageFile,folder,cond,reg,roundno)
{
	run("Bio-Formats Importer", "open=[" + folder + imageFile + "] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	//filename = getTitle();
	getDimensions(width, height, channels, slices, frames);

	Stack.setDisplayMode("composite"); //flattened version
	
	//run("Median...", "radius=3");
	run("Subtract Background...", "rolling=50");
	
	//UNCOMMENT THIS if you wish to run drift correction
	//run("Properties...", "channels="+ channels +"  slices=1 frames=" + (nSlices/4) + " pixel_width=0.1625000 pixel_height=0.1625000 voxel_depth=0.5000000 frame=4");
	//run("Correct 3D drift", "channel=1 correct only=0 lowest=1 highest=1 max_shift_x=50 max_shift_y=50 max_shift_z=10"); // drift correction, will determine if needed

	dotIndex = indexOf(imageFile, ".");	
	title = substring(imageFile, 0, dotIndex); 								
	titlesplit = split(title," ");
	roino = titlesplit[0];
	if (parseInt(roundno) < 10) {
		roundid = "round00" + roundno;
	} else if (parseInt(roundno) >= 10) {
		roundid = "round0" + roundno;
	}
	newname = cond + "-" + reg + "-" + roino + "_" + roundid;
	print(newname); //output file 

	saveAs("Tiff", parentout + newname + "_pp.tif");;

	close("*");
	run("Collect Garbage");
}
