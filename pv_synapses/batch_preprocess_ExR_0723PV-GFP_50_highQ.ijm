// ************* FIJI SCRIPT FOR PRE-PROCESSING EXR IMAGES ************* \\
// Last modified by MES on 8/22/24

run("Close All");
setBatchMode(true);

parent = "A:/Margaret/mExR/2024-08_SSC & HPC_PV-GFP + 5 Target_50step Crop Synapse ROI/"; //parent folder where raw microscope images are stored (e.g. in .nd2 format)
parentout = "A:/Margaret/mExR/2024-08_SSC & HPC_PV-GFP + 5 Target_50step Crop Synapse ROI/preprocessed/"; //where to save the images
//reg = "ctx"; //specify the brain region here (you can delete this from the final filename in the last section if you want)

folderList = newArray(
	"Cropped 50step stacks_Ctx/",
	"Cropped 50step stacks_HPC-CA1/"
	);
	
setList = newArray(
	"GluA1/",
	"GluA4/",
	"NR1/",
	"PSD95/",
	"Stg/"
	);
	
setBatchMode(true); //true for speed

for (m = 0; m < folderList.length; m++){
	for (j = 0; j < setList.length; j++){
	    	
		folder = parent+folderList[m]+setList[j];
		//print(folder);
		imgList = getFileList(folder);
		
	    	
		for (k = 0; k < imgList.length; k++){
	    	if (endsWith(imgList[k], ".tif")){ //CHANGE HERE for your filetype
	    		if (indexOf(imgList[k],"MAX") < 0){ //CHANGE HERE: provide a string that is ONLY in the images you're interested in
	    		
				//we have "40" because all of our high-mag images are named with "40x"
				// you can also skip this extra for-loop
					preprocessImage(imgList[k],folder);
	
				}
	    	}
		}
	}
}

function preprocessImage(imageFile,folder)
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


	saveAs("Tiff", parentout + title + "_pp.tif");;

	close("*");
	run("Collect Garbage");
}
