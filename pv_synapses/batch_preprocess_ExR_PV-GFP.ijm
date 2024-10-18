// ************* FIJI SCRIPT FOR PRE-PROCESSING MULTI-EXR IMAGES ************* \\
// Last modified by MES on 5/14/24
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

run("Close All");
setBatchMode(true);

parent = "A:/Jinyoung/2023.7.25 PV-GFP/"; //parent folder where raw microscope images are stored (e.g. in .nd2 format)
parentout = "A:/Margaret/mExR/2023.07_PV-GFP_ExR/"; //where to save the images
//reg = "ctx"; //specify the brain region here (you can delete this from the final filename in the last section if you want)

regList = newArray( 
	"CA1/",
	"SSC/"
	);

protList = newArray("GluA1/",
					"GluA4/",
					"NR1/",
					"PSD95/",
					"Stg/"); //round folders - again, will need to be modified depending on how you saved your data

sampleList = newArray(
	"A/",
	"B/",
	"C/"
	);
	
setBatchMode(true); //true for speed

for (m = 0; m < regList.length; m++){
	for (i = 0; i < protList.length; i++){
	    for (j = 0; j < sampleList.length; j++){
	    	
	    	folder = parent+regList[m]+protList[i]+sampleList[j];
	    	//print(folder);
	    	imgList = getFileList(folder);
	    	
	    	samplesplits = split(sampleList[j],"/");
	    	sample = samplesplits[0];

	    	protsplits = split(protList[i],"/");
	    	prot = protsplits[0];
	    	
	    	regsplits = split(regList[m],"/"); // uncomment this if you have region information in the directory structure, which will allow you to forgo hard-coding
			reg = regsplits[0];
	
			//CHANGE HERE per the naming convention of your files. This assumes "ROI[X] 40x.tif"
	    	//roundind = indexOf(roundList[j],"R");
	    	//slashind = indexOf(roundList[j],"/");
	    	//roundno = substring(roundList[j], roundind+1, slashind);
	    	roundno = toString(m + 1);
	    	
	    	for (k = 0; k < imgList.length; k++){
		    	if (endsWith(imgList[k], ".nd2")){ //CHANGE HERE for your filetype
		    		if (indexOf(imgList[k],"Max") < 0){ //CHANGE HERE: provide a string that is ONLY in the images you're interested in
		    			if (indexOf(imgList[k],"40x") >= 0){
		    			//we have "40" because all of our high-mag images are named with "40x"
		    			// you can also skip this extra for-loop
		    				preprocessImage(imgList[k],folder,sample,prot,reg);
		    			}
					}
				}
	    	}
		}
	}
}

function preprocessImage(imageFile,folder,sample,prot,reg)
{
	run("Bio-Formats Importer", "open=[" + folder + imageFile + "] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	//filename = getTitle();
	getDimensions(width, height, channels, slices, frames);

	//Stack.setDisplayMode("composite"); //flattened version
	
	//run("Median...", "radius=3");
	run("Subtract Background...", "rolling=50");
	
	//UNCOMMENT THIS if you wish to run drift correction
	//run("Properties...", "channels="+ channels +"  slices=1 frames=" + (nSlices/4) + " pixel_width=0.1625000 pixel_height=0.1625000 voxel_depth=0.5000000 frame=4");
	//run("Correct 3D drift", "channel=1 correct only=0 lowest=1 highest=1 max_shift_x=50 max_shift_y=50 max_shift_z=10"); // drift correction, will determine if needed

	dotIndex = indexOf(imageFile, ".");	
	title = substring(imageFile, 0, dotIndex); 								
	titlesplit = split(title,"x");
	roino = titlesplit[1];
	
	newname = reg + "_" + prot + "_" + sample + "_fov" + roino;
	print(newname); //output file 

	saveAs("Tiff", parentout + newname + "_pp.tif");;

	close("*");
	run("Collect Garbage");
}
