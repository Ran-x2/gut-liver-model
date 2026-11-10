# Prepare only a project-local R library and pinned official function sources.
args<-commandArgs(trailingOnly=TRUE)
script_arg<-grep('^--file=',commandArgs(),value=TRUE)
root<-if(length(args))normalizePath(args[[1]],winslash='/',mustWork=TRUE)else dirname(dirname(dirname(normalizePath(sub('^--file=','',script_arg[[1]]),winslash='/'))))
lib<-file.path(root,'runtime','R-library');dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=180)
cran<-c('dplyr','tibble','tidyr','ROCR','caTools','data.table','Matrix','ggplot2','digest','ragg')
missing<-cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,lib=lib,dependencies=NA,type=if(.Platform$OS.type=='windows')'binary'else'source')
if(!requireNamespace('limma',quietly=TRUE)){
 if(!requireNamespace('BiocManager',quietly=TRUE))install.packages('BiocManager',lib=lib)
 BiocManager::install('limma',lib=lib,ask=FALSE,update=FALSE,dependencies=NA)
}
required<-c(cran,'limma');stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))
commit<-'4ae630d1890894f37d921ab33c3760ebdde95b9c'
expected<-c('application_prediction.R'='d2012acb7697a35dc810465d28aec6b994a7fad587f7337e720b44fdf124763d','evaluate_model_ligand_prediction.R'='487f326221432385c6a3de60a3d9dc45afc02c35542fbccd9ee31e2b0d4447db','evaluate_model_target_prediction.R'='46c5515c02d4bb23f70df05e90696e37f8a0992b0ea12b792a4d3ea6f5e13420','supporting_functions.R'='5782390e932c0e0c89d38728a8178662a1837bfb9299d1538e666a8caf6abb5b')
vendor<-file.path(root,'work','nichenet_environment','nichenetr_v2.0.4');dir.create(vendor,recursive=TRUE,showWarnings=FALSE)
manifest<-list()
for(f in names(expected)){
 path<-file.path(vendor,f);url<-paste0('https://raw.githubusercontent.com/saeyslab/nichenetr/',commit,'/R/',f)
 if(!file.exists(path))download.file(url,path,mode='wb',quiet=TRUE)
 actual<-digest::digest(file=path,algo='sha256');if(actual!=expected[[f]])stop('Pinned source checksum mismatch: ',f)
 manifest[[f]]<-data.frame(file=f,tag='v2.0.4',commit=commit,url=url,sha256=actual,license='GPL-3; official unmodified source')
}
write.table(do.call(rbind,manifest),file.path(root,'work','nichenet_environment','official_source_manifest.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
info<-data.frame(package=required,version=vapply(required,function(x)as.character(packageVersion(x)),character(1)),path=vapply(required,find.package,character(1)))
write.table(info,file.path(root,'work','nichenet_environment','dependency_versions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
source(file.path(root,'scripts','nichenet','validate_official_functions.R'),local=FALSE)
cat('Prepared project-local dependencies and verified official source checksums.\n')