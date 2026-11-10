# Pinned third-party source retains GPL-3; see vendor/LICENSE-GPL-3.txt.
# Reuse recovered DE; reconstruct input/background tables; correct only NicheNet scoring.
# Exact NicheNet v2.0.4 functions are sourced unmodified from the pinned official commit.
args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
root <- if(length(args)) normalizePath(args[[1]],winslash='/',mustWork=TRUE) else dirname(dirname(dirname(normalizePath(sub('^--file=','',script_arg[[1]]),winslash='/'))))
.libPaths(c(Sys.getenv('NICHE_R_LIBRARY'),file.path(root,'runtime','R-library'),.libPaths()))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
options(stringsAsFactors=FALSE)
set.seed(20261110)
output_root<-Sys.getenv('REVISION_OUTPUT_DIR',unset=file.path(root,'output'))
out6<-file.path(output_root,'06_nichenet_inputs');out7<-file.path(output_root,'07_nichenet_activity');out8<-file.path(output_root,'08_receiver_targets')
out6_results<-file.path(out6,'output');out6_audit<-file.path(out6,'audit')
out7_results<-file.path(out7,'output');out7_audit<-file.path(out7,'audit')
out8_results<-file.path(out8,'output')
export_score_comparisons<-tolower(Sys.getenv('NICHE_EXPORT_SCORE_COMPARISONS',unset='false'))%in%c('1','true','yes')
export_input_details<-tolower(Sys.getenv('NICHE_EXPORT_INPUT_DETAILS',unset='false'))%in%c('1','true','yes')
for(p in c(out6,out6_results,out6_audit,out7,out7_results,out7_audit,out8,out8_results,file.path(out6,'gene_lists'),file.path(out6,'backgrounds')))dir.create(p,recursive=TRUE,showWarnings=FALSE)
write_tsv<-function(x,path)write.table(x,path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
sha<-function(path)digest::digest(file=path,algo='sha256')
log_path<-file.path(out7_audit,'nichenet_activity_execution_log.txt');writeLines(character(),log_path)
log_note<-function(...) {line<-paste(...,collapse=' ');cat(line,'\n');cat(line,'\n',file=log_path,append=TRUE)}
deg_dir<-file.path(root,'data','original_deg');cache_dir<-file.path(root,'data','original_nichenet')
model_path<-Sys.getenv('NICHE_MODEL');if(!nzchar(model_path))stop('Set NICHE_MODEL.')
lr_path<-Sys.getenv('NICHE_LR');if(!nzchar(lr_path))stop('Set NICHE_LR.')
source_files<-c('supporting_functions.R','evaluate_model_target_prediction.R','evaluate_model_ligand_prediction.R','application_prediction.R')
source_commit<-'4ae630d1890894f37d921ab33c3760ebdde95b9c'
vendor<-Sys.getenv('NICHE_VENDOR',unset=file.path(root,'work','nichenet_environment','nichenetr_v2.0.4'))
if(!all(file.exists(file.path(vendor,source_files))))stop('Prepare the pinned official NicheNet functions before execution.')
nn<-new.env(parent=globalenv());for(f in source_files)sys.source(file.path(vendor,f),envir=nn)
source_manifest<-data.frame(file=source_files,commit=source_commit,sha256=vapply(file.path(vendor,source_files),sha,character(1)),url=paste0('https://raw.githubusercontent.com/saeyslab/nichenetr/',source_commit,'/R/',source_files))
write_tsv(source_manifest,file.path(out7_audit,'nichenet_function_source_manifest.tsv'))
model_sha<-sha(model_path);lr_sha<-sha(lr_path)
if(model_sha!='699fce17ff65e2511277359696306bb130aa833bd7f30dd9290ad9efd9dc9c5d')stop('Candidate model differs from the audited reference; re-audit model matching.')
mat<-readRDS(model_path);lr<-readRDS(lr_path)
stopifnot(is.matrix(mat),!anyDuplicated(rownames(mat)),all(c('from','to')%in%names(lr)))
model_genes<-rownames(mat)
write_tsv(data.frame(role=c('candidate_ligand_target_matrix','candidate_ligand_receptor_network'),path=c(model_path,lr_path),sha256=c(model_sha,lr_sha),note=c('Recovered candidate reference; original version identity not assumed.','Recovered candidate reference; original version identity not assumed.')),file.path(out6_audit,'reference_nichenet_manifest.tsv'))
receivers<-c('gutT','gutMac','liverT','liverMac')
receiver_labels<-c(gutT='Gut T cells',gutMac='Gut macrophages',liverT='Liver T cells',liverMac='Liver macrophages')
receiver_file_labels<-c(gutT='gut_T_cells',gutMac='gut_macrophages',liverT='liver_T_cells',liverMac='liver_macrophages')
donor_names<-c('3'='Donor AJD3280','4'='Donor AJG2309');donor_labels<-c('3'='Donor 1','4'='Donor 2')
expected_post<-c('3_gutT'=484L,'4_gutT'=162L,'3_gutMac'=6L,'4_gutMac'=42L,'3_liverT'=1133L,'4_liverT'=823L,'3_liverMac'=1068L,'4_liverMac'=1232L)
degs<-list();backgrounds<-list();gene_sets<-list();case_meta<-list();input_manifest<-list();schema<-list();background_checks<-list()
for(receiver in receivers){
  csv_path<-file.path(deg_dir,paste0(receiver,'_DEG.csv'));obj_path<-file.path(deg_dir,paste0(receiver,'.gz'))
  if(!all(file.exists(c(csv_path,obj_path))))stop('Missing recovered original inputs for ',receiver)
  dat<-read.table(csv_path,sep=',',header=TRUE,row.names=1,stringsAsFactors=FALSE)
  stopifnot(all(c('gene_short_name','term','q_value','normalized_effect')%in%names(dat)),all(dat$term=='DonorIDDonor AJG2309'),!anyNA(dat),!anyDuplicated(dat$gene_short_name),all(nzchar(dat$gene_short_name)),all(is.finite(dat$q_value)),all(dat$q_value>=0&dat$q_value<=1),all(is.finite(dat$normalized_effect)))
  degs[[receiver]]<-dat
  obj<-readRDS(obj_path);md<-attr(obj,'meta.data');assay<-attr(obj,'assays')[['RNA']];cm<-attr(assay,'layers')[['counts']]
  genes<-attr(attr(assay,'features'),'dimnames')[[1]];cells<-attr(attr(assay,'cells'),'dimnames')[[1]]
  stopifnot(inherits(cm,'dgCMatrix'),identical(dim(cm),c(length(genes),length(cells))),identical(genes,dat$gene_short_name),!anyDuplicated(cells),all(cells%in%rownames(md)))
  md<-md[match(cells,rownames(md)),,drop=FALSE];stopifnot(identical(cells,rownames(md)),is.factor(md$DonorID),identical(levels(md$DonorID),unname(donor_names)))
  stopifnot(all(is.finite(cm@x)),all(cm@x>=0))
  rownames(cm)<-genes;colnames(cm)<-cells
  schema[[receiver]]<-data.frame(receiver=receiver,table_rows=nrow(dat),unique_genes=length(unique(genes)),receiver_cells=ncol(cm),donor_factor_levels=paste(levels(md$DonorID),collapse=';'),invalid_values=0L,duplicate_genes=0L)
  input_manifest[[paste0(receiver,'_csv')]]<-data.frame(receiver=receiver,role='original_saved_DEG_table',path=csv_path,bytes=file.info(csv_path)$size,sha256=sha(csv_path))
  input_manifest[[paste0(receiver,'_rds')]]<-data.frame(receiver=receiver,role='original_DE_input_Seurat_object',path=obj_path,bytes=file.info(obj_path)$size,sha256=sha(obj_path))
  for(donor_code in c('3','4')){
    case_id<-paste(donor_code,receiver,sep='_');cols<-which(as.character(md$DonorID)==donor_names[[donor_code]]);n_cells<-length(cols);stopifnot(n_cells>0)
    donor_counts<-cm[,cols,drop=FALSE]
    n_expr<-Matrix::rowSums(donor_counts>0);frac<-n_expr/n_cells;bg<-genes[frac>=0.01];backgrounds[[case_id]]<-bg
    positive<-if(donor_code=='3')dat$q_value<=0.05&dat$normalized_effect<=-0.5 else dat$q_value<=0.05&dat$normalized_effect>=0.5
    pre<-dat[positive,,drop=FALSE];matched<-pre$gene_short_name%in%model_genes;post<-pre[matched,,drop=FALSE];gene_sets[[case_id]]<-post$gene_short_name
    stopifnot(nrow(post)==expected_post[[case_id]])
    post$directional_effect<-if(donor_code=='3')-post$normalized_effect else post$normalized_effect
    write_tsv(post,file.path(out6,'gene_lists',paste0(case_id,'_in_nichenet.tsv')))
    if(export_input_details){
      write_tsv(pre,file.path(out6,'gene_lists',paste0(case_id,'_pre_nichenet.tsv')))
      write_tsv(pre[!matched,,drop=FALSE],file.path(out6,'gene_lists',paste0(case_id,'_not_in_nichenet.tsv')))
    }
    write_tsv(data.frame(gene_short_name=bg),file.path(out6,'backgrounds',paste0(case_id,'_expressed_background.tsv')))
    if(export_input_details){
      detected<-data.frame(gene_short_name=genes,n_cells_expressing=n_expr,n_cells=n_cells,fraction_expressing=frac,mean_source_count=as.numeric(Matrix::rowMeans(donor_counts)))
      write_tsv(detected,file.path(out6,'backgrounds',paste0(case_id,'_expression_detection.tsv')))
    }
    oldbg_path<-file.path(cache_dir,paste0(case_id,'_receptors.csv'))
    if(file.exists(oldbg_path)){
      oldbg<-read.csv(oldbg_path,row.names=1,check.names=FALSE)[[1]]
      background_checks[[case_id]]<-data.frame(case_id=case_id,cached_background_genes=length(oldbg),reconstructed_background_genes=length(bg),identical_order=identical(oldbg,bg),set_equal=setequal(oldbg,bg),only_in_cache=length(setdiff(oldbg,bg)),only_in_reconstruction=length(setdiff(bg,oldbg)))
      if(!setequal(oldbg,bg))stop('Cached original background disagrees with recovered DE input object: ',case_id)
    }
    bg_model<-intersect(bg,model_genes);union_genes<-union(bg_model,post$gene_short_name)
    case_meta[[case_id]]<-data.frame(case_id=case_id,receiver=receiver,receiver_label=receiver_labels[[receiver]],donor_code=donor_code,donor_label=donor_labels[[donor_code]],donor_name=donor_names[[donor_code]],receiver_cells=n_cells,effect_cutoff=if(donor_code=='3')'<=-0.5' else '>=0.5',q_cutoff='<=0.05',DEGs_before_model=nrow(pre),DEGs_after_model=nrow(post),unique_DEGs_after_model=length(unique(post$gene_short_name)),removed_by_model=sum(!matched),expressed_background_genes=length(bg),background_model_genes=length(bg_model),positive_genes_outside_background=length(setdiff(post$gene_short_name,bg_model)),scored_gene_universe=length(union_genes),positive_prevalence=nrow(post)/length(union_genes),low_information=nrow(post)<20L)
  }
  rm(obj,assay,cm);invisible(gc())
}
meta<-do.call(rbind,case_meta);rownames(meta)<-NULL
write_tsv(meta,file.path(out6_results,'nichenet_input_counts.tsv'));write_tsv(do.call(rbind,schema),file.path(out6_audit,'input_schema_audit.tsv'));write_tsv(do.call(rbind,input_manifest),file.path(out6_audit,'original_input_manifest.tsv'))
if(length(background_checks))write_tsv(do.call(rbind,background_checks),file.path(out6_audit,'cached_background_concordance.tsv'))
log_note('Original DE tables and input objects audited. No DE, normalization or annotation rerun.')
activities<-list();cache_validation<-list();correction_comparisons<-list();target_links<-list();display_targets<-list();ligand_support<-list()
for(case_id in meta$case_id){
  info<-meta[meta$case_id==case_id,,drop=FALSE];receiver<-info$receiver;donor_code<-info$donor_code;bg<-backgrounds[[case_id]];gs<-gene_sets[[case_id]]
  potential<-unique(lr$from[lr$to%in%bg]);stopifnot(length(potential)>0,all(potential%in%colnames(mat)))
  cache_path<-file.path(cache_dir,paste0(case_id,'_nichenet.csv'));cached<-NULL;reuse<-FALSE
  if(file.exists(cache_path)){
    cached<-read.table(cache_path,header=TRUE,stringsAsFactors=FALSE)
    stopifnot(all(c('test_ligand','auroc','aupr','aupr_corrected','pearson')%in%names(cached)),!anyDuplicated(cached$test_ligand))
    correct_candidate_match<-setequal(cached$test_ligand,potential)
    correct_prevalence_error<-max(abs((cached$aupr-cached$aupr_corrected)-info$positive_prevalence))
    legacy_bg<-backgrounds[[paste('3',receiver,sep='_')]];legacy_potential<-unique(lr$from[lr$to%in%legacy_bg]);legacy_prev<-length(gs)/length(union(intersect(legacy_bg,model_genes),gs))
    legacy_prevalence_error<-max(abs((cached$aupr-cached$aupr_corrected)-legacy_prev))
    probe_error<-NA_real_
    if(donor_code=='3'&&correct_candidate_match&&correct_prevalence_error<1e-10){
      probe_ligands<-unique(c(cached$test_ligand[[1]],cached$test_ligand[[ceiling(nrow(cached)/2)]],tail(cached$test_ligand,1)))
      probe<-nn$predict_ligand_activities(gs,bg,mat,probe_ligands)
      cache_probe<-cached[match(probe$test_ligand,cached$test_ligand),]
      probe_values<-as.matrix(probe[,c('auroc','aupr','aupr_corrected','pearson')]);cache_values<-as.matrix(cache_probe[,c('auroc','aupr','aupr_corrected','pearson')]);probe_error<-if(!identical(is.na(probe_values),is.na(cache_values)))Inf else max(c(0,abs(probe_values-cache_values)),na.rm=TRUE)
      reuse<-is.finite(probe_error)&&probe_error<1e-9
    }
    cache_validation[[case_id]]<-data.frame(case_id=case_id,cached_ligands=nrow(cached),corrected_candidate_set_matches=correct_candidate_match,corrected_prevalence_error=correct_prevalence_error,legacy_donor3_candidate_set_matches=setequal(cached$test_ligand,legacy_potential),legacy_donor3_prevalence_error=legacy_prevalence_error,three_ligand_probe_max_error=probe_error,reused_original_donor3_scores=reuse,cached_sha256=sha(cache_path))
  }
  if(reuse){act<-cached;mode<-'recovered_original_donor3_scores_validated_on_three_probes';log_note(case_id,': reused cached donor3 table after candidate, prevalence and three-ligand validation.')}else{log_note(case_id,': evaluating',length(potential),'ligands with corrected donor-specific background.');act<-as.data.frame(nn$predict_ligand_activities(gs,bg,mat,potential));mode<-if(donor_code=='4')'rescored_with_correct_donor4_background' else 'scored_missing_original_donor3_table'}
  act<-act[order(-act$auroc,act$test_ligand),];act$rank_auroc<-seq_len(nrow(act));act$rank_aupr_corrected<-rank(-act$aupr_corrected,ties.method='min');act$case_id<-case_id;act$receiver<-receiver;act$donor_code<-donor_code;act$score_origin<-mode
  stopifnot(all(is.finite(act$auroc)),all(act$auroc>=0&act$auroc<=1))
  activities[[case_id]]<-act
  if(!is.null(cached)&&donor_code=='4'){
    comparison_path<-file.path(out7_audit,paste0('donor2_',receiver_file_labels[[receiver]],'_original_vs_corrected_ligand_scores.tsv'))
    if(export_score_comparisons){
      cmp<-merge(cached,act,by='test_ligand',all=TRUE,suffixes=c('_original','_corrected'))
      cmp$case_id<-case_id;cmp$auroc_change<-cmp$auroc_corrected-cmp$auroc_original;write_tsv(cmp,comparison_path)
    }
    common<-intersect(cached$test_ligand,act$test_ligand)
    correction_comparisons[[case_id]]<-data.frame(case_id=case_id,original_candidates=nrow(cached),corrected_candidates=nrow(act),gained_candidates=length(setdiff(act$test_ligand,cached$test_ligand)),lost_candidates=length(setdiff(cached$test_ligand,act$test_ligand)),max_absolute_AUROC_change=max(abs(cached$auroc[match(common,cached$test_ligand)]-act$auroc[match(common,act$test_ligand)])),original_top_ligand=cached$test_ligand[[which.max(cached$auroc)]],corrected_top_ligand=act$test_ligand[[1]])
  }
  top_ligands<-head(act$test_ligand,15);top3<-head(top_ligands,3)
  de<-degs[[receiver]];de<-de[match(gs,de$gene_short_name),];de$directional_effect<-if(donor_code=='3')-de$normalized_effect else de$normalized_effect
  for(ligand in top_ligands){
    links<-as.data.frame(nn$get_weighted_ligand_target_links(ligand,gs,mat,n=250))
    links<-links[!is.na(links$target),,drop=FALSE]
    ligand_support[[paste(case_id,ligand)]]<-data.frame(case_id=case_id,ligand=ligand,AUROC=act$auroc[match(ligand,act$test_ligand)],receiver_DEGs_among_top250_prior_targets=nrow(links),target_prior_cutoff=250L,input_DEGs=length(gs))
    if(nrow(links)){
      links<-cbind(data.frame(case_id=case_id,receiver=receiver,donor_code=donor_code),links,de[match(links$target,de$gene_short_name),c('q_value','normalized_effect','directional_effect'),drop=FALSE])
      links<-links[order(-links$weight,links$target),];target_links[[paste(case_id,ligand)]]<-links
      if(ligand%in%top3)display_targets[[paste(case_id,ligand)]]<-head(links,5)
    }
  }
}
all_activity<-do.call(rbind,activities);rownames(all_activity)<-NULL
write_tsv(all_activity,file.path(out7_results,'ligand_activity_rankings_all_donors_and_cell_types.tsv'))
if(length(cache_validation))write_tsv(do.call(rbind,cache_validation),file.path(out7_audit,'cached_ligand_score_validation.tsv'))
if(length(correction_comparisons))write_tsv(do.call(rbind,correction_comparisons),file.path(out7_audit,'donor2_background_correction_summary.tsv'))
all_links<-do.call(rbind,target_links);rownames(all_links)<-NULL
write_tsv(all_links,file.path(out8_results,'upregulated_receiver_gene_links_for_top15_ligands.tsv'));write_tsv(do.call(rbind,ligand_support),file.path(out8_results,'upregulated_receiver_gene_counts_for_top15_ligands.tsv'))
display<-do.call(rbind,display_targets);rownames(display)<-NULL
# Plotting subsets stay in memory; complete links and counts are final deliverables.
# Shared figure style: Arial, white background, 600-dpi PNG only.
theme_revision<-function()theme_classic(base_family='Arial',base_size=10)+theme(plot.title=element_text(size=14,face='bold',hjust=0.5,margin=margin(b=8)),plot.title.position='plot',plot.subtitle=element_blank(),axis.title=element_text(size=10),axis.text=element_text(size=8,color='#222222'),strip.background=element_rect(fill='#F2F4F6',color='#CCD3D9',linewidth=0.4),strip.text=element_text(size=9,face='bold',lineheight=1.1,margin=margin(5,5,5,5)),legend.position='bottom',legend.title=element_text(size=9),legend.text=element_text(size=8),panel.spacing=grid::unit(10,'pt'),plot.margin=margin(12,15,12,12))
save_figure<-function(plot,path,width,height){
  ragg::agg_png(paste0(path,'.png'),width=width,height=height,units='in',res=600,background='white');print(plot);grDevices::dev.off()
}
meta$panel_label<-paste0(meta$receiver_label,' | ',meta$donor_label,'\n',format(meta$DEGs_after_model,big.mark=',',trim=TRUE),' input DEGs; ',format(meta$receiver_cells,big.mark=',',trim=TRUE),' cells')
panel_levels<-meta$panel_label
add_panel<-function(df){df$panel_label<-factor(meta$panel_label[match(df$case_id,meta$case_id)],levels=panel_levels);df$donor_label<-factor(donor_labels[df$donor_code],levels=c('Donor 1','Donor 2'));df}
all_activity<-add_panel(all_activity)
top<-all_activity[all_activity$rank_auroc<=15,,drop=FALSE]
top$row_id<-factor(paste(top$case_id,top$test_ligand,sep='___'),levels=rev(paste(top$case_id,top$test_ligand,sep='___')))
donor_colors<-c('Donor 1'='#3B6484','Donor 2'='#A46A51')
p_activity<-ggplot(top,aes(x=auroc,y=row_id,color=donor_label))+geom_vline(xintercept=0.5,linetype='dashed',color='#A7ADB4',linewidth=0.4)+geom_point(size=2.1)+geom_text(aes(label=sprintf('%.4f',auroc)),hjust=-0.25,size=2.45,family='Arial',show.legend=FALSE)+facet_wrap(~panel_label,ncol=2,scales='free_y')+scale_y_discrete(labels=function(x)sub('^.*___','',x))+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),expand=expansion(mult=c(0.01,0.02)))+scale_color_manual(values=donor_colors,guide='none')+labs(title='Candidate ligands prioritized by receiver DEGs',x='NicheNet ligand activity (AUROC)',y=NULL)+theme_revision()
save_figure(p_activity,file.path(out7_results,'Figure_R7A_top15_ligands_shared_AUROC_scale'),10.5,12.5)
p_distribution<-ggplot(all_activity,aes(x=auroc,color=donor_label))+stat_ecdf(geom='step',linewidth=0.7)+geom_vline(xintercept=0.5,linetype='dashed',color='#A7ADB4',linewidth=0.4)+facet_wrap(~panel_label,ncol=2)+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),expand=expansion(mult=0))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.25),expand=expansion(mult=c(0,0.02)))+scale_color_manual(values=donor_colors,guide='none')+labs(title='Distribution of all candidate-ligand scores',x='NicheNet ligand activity (AUROC)',y='Cumulative fraction of ligands')+theme_revision()+theme(panel.spacing.x=grid::unit(26,'pt'))
save_figure(p_distribution,file.path(out7_results,'Figure_R7B_all_ligand_AUROC_distributions'),10.5,9.5)
display<-add_panel(display)
display$pair_label<-paste(display$ligand,'\u2192',display$target)
display$row_key<-paste(display$case_id,display$pair_label,sep='___')
display$q_category<-factor(ifelse(display$q_value<=0.001,'q <= 0.001',ifelse(display$q_value<=0.01,'q <= 0.01','q <= 0.05')),levels=c('q <= 0.001','q <= 0.01','q <= 0.05'))
missing_cases<-meta[!meta$case_id%in%display$case_id,,drop=FALSE]
plot_targets<-display
if(nrow(missing_cases)){
  dummy<-display[rep(1L,nrow(missing_cases)),,drop=FALSE]
  dummy$case_id<-missing_cases$case_id;dummy$panel_label<-factor(missing_cases$panel_label,levels=panel_levels);dummy$directional_effect<-NA_real_;dummy$weight<-NA_real_;dummy$row_key<-paste(missing_cases$case_id,'No supported targets',sep='___');dummy$pair_label<-'No supported targets'
  plot_targets<-rbind(display,dummy)
}
plot_targets$row_key<-factor(plot_targets$row_key,levels=rev(unique(plot_targets$row_key)))
x_upper<-ceiling(max(display$directional_effect)*1.1*2)/2
p_targets<-ggplot(plot_targets,aes(x=directional_effect,y=row_key))+geom_vline(xintercept=0.5,linetype='dashed',color='#A7ADB4',linewidth=0.4)+geom_point(aes(size=weight,color=q_category),alpha=0.9,na.rm=TRUE)+facet_wrap(~panel_label,ncol=2,scales='free_y',drop=FALSE)+scale_y_discrete(labels=function(x)sub('^.*___','',x))+scale_x_continuous(limits=c(0,x_upper),expand=expansion(mult=c(0.01,0.03)))+scale_color_manual(values=c('q <= 0.001'='#204F70','q <= 0.01'='#6E95B0','q <= 0.05'='#B7CBD9'),drop=FALSE,name='Saved receiver DE q value')+scale_size_continuous(range=c(1.8,4.1),name='Model regulatory potential')+labs(title='Upregulated receiver genes linked to prioritized ligands',x='Monocle3 normalized_effect toward indicated donor',y='Ligand \u2192 receiver target')+theme_revision()+theme(legend.position='bottom',legend.box='vertical')
if(nrow(missing_cases))p_targets<-p_targets+geom_text(data=plot_targets[!is.finite(plot_targets$directional_effect),,drop=FALSE],aes(x=x_upper*0.52,y=row_key),label='No input DEGs among these\nligands\u2019 top 250 model targets',inherit.aes=FALSE,size=3,family='Arial',color='#454545')
save_figure(p_targets,file.path(out8_results,'Figure_R8_ligand_links_to_upregulated_receiver_genes'),11.5,13)
range_rows<-lapply(meta$case_id,function(id){a<-activities[[id]];t<-head(a$auroc,15);data.frame(case_id=id,candidate_ligands=nrow(a),top_ligand=a$test_ligand[[1]],top_AUROC=t[[1]],top15_min=min(t),top15_max=max(t),top15_range=max(t)-min(t),score_origin=a$score_origin[[1]])})
ranges<-do.call(rbind,range_rows);write_tsv(ranges,file.path(out7_results,'top_ligands_and_score_ranges_by_donor_and_cell_type.tsv'))
writeLines(c("Figure: Candidate ligands prioritized by receiver DEGs. Each panel shows the 15 highest AUROC-scored candidate ligands for one receiver population and donor direction. All panels use the same 0-1 AUROC axis; the dashed line marks AUROC 0.5. Printed values retain four decimal places. Candidate ligands were selected by receptor expression in the receiver, with expression defined as counts >0 in at least 1% of the original donor-specific receiver cells. The original receiver DEG tables were reused; q_value <=0.05 and directional normalized_effect of at least0.5 were applied before matching the model gene universe. Donor 2 analyses use their own corrected expressed background instead of the Donor 1 background used by the deposited script. Scores rank candidate upstream ligands and are not evidence of uniquely causal signaling or of ligand secretion. Donor 1 gut macrophages have only six input DEGs (RPS8, EEF1A1, RPS26, RPS2, RPL13, RPS4X); interpret that panel descriptively.","Complete scores and rankings for all candidate ligands are provided in the accompanying table."),file.path(out7_results,'Figure_R7A_top15_ligands_shared_AUROC_scale_legend.txt'))
writeLines(c('Figure: Distribution of all candidate-ligand scores. Empirical cumulative distribution curves show the fraction of candidate ligands with AUROC at or below each value, for all eight donor/receiver cases. Both axes are fixed0-1 across panels. The dashed line marks AUROC0.5. This descriptive display contextualizes the narrow ranges among the highest-ranked ligands; no hypothesis test comparing ligands was performed. The input DEG counts and receiver cell counts are shown in each panel.'),file.path(out7_results,'Figure_R7B_all_ligand_AUROC_distributions_legend.txt'))
writeLines(c("Figure to reviewer comment 8. Upregulated receiver genes linked to prioritized ligands.","For the three highest-AUROC ligands in each donor/cell-type case, the receiver DEG set was intersected with each ligand's top 250 NicheNet prior targets. Up to five matching links per ligand, ordered by regulatory potential, are displayed. The accompanying tables provide complete matching links and receiver-gene counts, including zero-support cases, for the top 15 ligands.","The x-axis shows the saved Monocle3 normalized_effect oriented toward the indicated donor; the dashed line marks the 0.5 input threshold. Point size represents NicheNet prior regulatory potential, and color represents the saved receiver-gene q-value category. A panel with no matching input DEGs explicitly indicates the absence of supported targets under these criteria.","The receiver-gene evidence comes from the same DEGs used to prioritize the ligands. It therefore describes transcriptional consistency rather than independent validation of ligand expression, secretion, or causal signaling."),file.path(out8_results,'Figure_R8_ligand_links_to_upregulated_receiver_genes_legend.txt'))
writeLines(c('Implementation: official NicheNet functions from tagv2.0.4, pinned commit4ae630d1890894f37d921ab33c3760ebdde95b9c, sourced unmodified. This avoids installing unrelated random-forest, optimization and Seurat dependencies. Metric dependencies are installed only in the staging-local R library. Six synthetic validation checks, including an independent Mann-Whitney AUROC calculation, are recorded in work/nichenet_environment.', '', 'Changes to deposited workflow: donor4 uses its own expressed background; candidate scores are sorted explicitly; active ligand-target links are restricted to observed receiver DEGs; plot scales are fixed. Original saved receiver DE tables and original receiver cohorts are retained. No DE fitting, normalization, clustering or annotation is rerun.', '', 'Original donor3 cached activity tables are reused only after candidate-set and positive-prevalence concordance plus three independently rescored ligand probes agree within1e-9. Missing donor3 tables and donor4 corrections are scored with the official functions. Every row records score_origin. Target-link display uses top15 ligands by AUROC, the standard top250 prior-target intersection, and up tofive targets for each topthree ligands for visualization.', '', 'The manuscript section lists NicheNetv2.2.0 whereas the deposited script comments specifyv2.0.4. This implementation records its exact source and validates against recovered scores; it does not claim that an undocumented historical package installation was recovered.'),file.path(out7_audit,'nichenet_software_and_analysis_scope.txt'))
capture.output(sessionInfo(),file=file.path(out7_audit,'R_session_information.txt'))
log_note('Completed input audit, corrected ligand activity and receiver target support. No DE refits were performed.')
