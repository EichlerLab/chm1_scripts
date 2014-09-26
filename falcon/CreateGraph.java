



//Simple PDF export
ExportController ec = Lookup.getDefault().lookup(ExportController.class);
try {
   ec.exportFile(new File("simple.pdf"));
} catch (IOException ex) {
   ex.printStackTrace();
   return;
}
 
//PDF Exporter config and export to Byte array
PDFExporter pdfExporter = (PDFExporter) ec.getExporter("pdf");
pdfExporter.setPageSize(PageSize.A0);
ByteArrayOutputStream baos = new ByteArrayOutputStream();
ec.exportStream(baos, pdfExporter);
byte[] pdf = baos.toByteArray();
