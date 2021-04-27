// input files location
var inputFolder = Folder.selectDialog("Choose the folder that contains psd files");
var inputFiles = inputFolder.getFiles("*.psd");

// output files location
var outputFolder = Folder.selectDialog("Choose the folder to save files to");

// SET the options for saving pdf files
var PDFOptions = new PDFSaveOptions();
PDFOptions.PDFCompatibility = PDFCompatibility.PDF13


for(var i = 0; i<inputFiles.length; i++){
	open(inputFiles[i]);
	app.activeDocument.saveAs(
	new File(outputFolder + "/" + inputFiles[i].name + ".pdf"),
	PDFOptions
	);
	app.activeDocument.close();
	// alert(inputFiles[i].name)
}
alert("All psd Files has been converted to PDF!");