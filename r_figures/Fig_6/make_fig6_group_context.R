file_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]),winslash="/",mustWork=TRUE))
repo <- normalizePath(file.path(script_dir,"../.."),winslash="/",mustWork=TRUE)
ROOT <- file.path(repo,'mechanism')
source(file.path(repo,'r_figures/common/mechanism_style.R'))
D <- Sys.getenv('DMWA_W_DATA_DIR',unset=file.path(ROOT,'data','w'))
rd <- function(name) read_csv(file.path(D,name),show_col_types=FALSE)

ini <- rd('toy_initial_states.csv') |> filter(network=='toy_1') |>
  mutate(group=factor(members,levels=c('AB','AC','BC','AD','AE','ABC')),
         state=factor(state,levels=c(1,2),labels=c('Triangle state','Edge states')),
         xpos=as.numeric(group)+ifelse(state=='Triangle state',-.12,.12))
pA <- ggplot(ini,aes(xpos,mass,colour=state,shape=state))+
  geom_segment(aes(xend=xpos,y=0,yend=mass),linewidth=.65)+
  geom_point(fill='white',size=2.4,stroke=.8)+
  scale_colour_manual(values=c('Triangle state'=BLUE,'Edge states'=ORANGE))+
  scale_shape_manual(values=c('Triangle state'=21,'Edge states'=24))+
  scale_x_continuous(breaks=1:6,labels=levels(ini$group),expand=expansion(mult=c(.06,.06)))+
  scale_y_continuous(breaks=c(0,1/3,2/3,1),labels=c('0','1/3','2/3','1'),expand=expansion(mult=c(.02,.06)))+
  annotate('text',x=2.65,y=.88,label='C*z[1]==C*z[2]',parse=TRUE,family='Arial',size=FIG_GEOM_TEXT_SIZE)+
  labs(title=ttl('A','Same nodes, different group states'),x='Hyperedge',y='Initial group mass')+
  theme(panel.grid.major.x=element_blank(),legend.position='bottom',legend.text=element_text(size=9.8),legend.key.width=unit(6,'mm'))

cl <- rd('closure_summary.csv') |>
  mutate(variant=factor(variant,levels=c('Full DMWA','Constant bias','Coverage'),labels=c('Full W','Constant bias','Coverage-only')),
         context=factor(projection,levels=c('Node','Node + order','Node + order + degree'),labels=c('Nodes','Nodes +\norder','Nodes +\norder + degree')),
         value=ifelse(maximum<1e-10,0,median),label=ifelse(maximum<1e-10,'Exact',sprintf('%.3f',median)))
pB <- ggplot(cl,aes(context,variant,fill=value))+
  geom_tile(width=.96,height=.93,colour='white',linewidth=.6)+
  geom_text(aes(label=label,colour=value>.06),family='Arial',size=FIG_GEOM_TEXT_SIZE)+
  scale_fill_gradient(low='#F0F3F7',high=BLUE,limits=c(0,max(cl$value)),guide='none')+
  scale_colour_manual(values=c('FALSE'=INK,'TRUE'='white'),guide='none')+
  labs(title=ttl('B','Closure by retained context'),x=NULL,y=NULL)+
  theme(panel.grid=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(size=9.8))

sep <- rd('same_node_separation.csv') |> filter(!is_toy,variant=='Full DMWA')
pC <- ggplot(sep,aes(100*total_fraction,spatial_TV_mean))+
  geom_point(shape=21,fill='white',colour=BLUE,size=1.9,stroke=.6,alpha=.8)+
  scale_x_continuous(breaks=c(0,1,2),expand=expansion(mult=c(.04,.08)))+
  scale_y_continuous(labels=function(x)sprintf('%.3f',x),n.breaks=4,expand=expansion(mult=c(.04,.08)))+
  labs(title=ttl('C','Total and spatial response separation'),x='Total response separation (%)',y='Mean spatial distance (TV)')

err <- rd('response_errors.csv') |> filter(!is_toy,variant=='Full DMWA') |>
  mutate(context=factor(projection,levels=c('Node','Node + order','Node + order + degree'),labels=c('Nodes','Nodes +\norder','Nodes +\norder + degree')),
         error_percent=100*relative_L2_error)
es <- err |> group_by(context) |> summarise(median=median(error_percent),q25=quantile(error_percent,.25),q75=quantile(error_percent,.75),.groups='drop') |>
  mutate(label=ifelse(median<1e-10,'Numerical\nzero',sprintf('%.3f%%',median)))
pD <- ggplot(err,aes(context,error_percent,colour=context))+
  geom_point(position=position_jitter(width=.12,height=0,seed=260916),shape=21,fill='white',size=1.3,stroke=.4,alpha=.5)+
  geom_errorbar(data=es,aes(y=median,ymin=q25,ymax=q75),width=.14,linewidth=.7)+
  geom_point(data=es,aes(y=median),shape=21,fill='white',size=2.5,stroke=.9)+
  geom_text(data=es,aes(y=7.1,label=label),family='Arial',size=FIG_GEOM_TEXT_SIZE,colour=INK,lineheight=.95)+
  scale_colour_manual(values=c('Nodes'=BLUE,'Nodes +\norder'=ORANGE,'Nodes +\norder + degree'=GREEN),guide='none')+
  scale_y_continuous(breaks=c(0,2,4,6),limits=c(-.15,7.8),expand=expansion(mult=0))+
  labs(title=ttl('D','Recovering node-response trajectories'),x=NULL,y='Relative response error (%)')+
  theme(axis.text.x=element_text(size=9.8),panel.grid.major.x=element_blank())

# Local layout revision only; source values and panel encodings are unchanged.
OUT <- Sys.getenv('DMWA_FIG6_OUT',unset=file.path(repo,'outputs','Fig_6'))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
p <- wrap_plots(pA,pB,pC,free(pD,side='l'),ncol=2)
# Figure contract: group order and overlap degree recover projected dynamics.
# A: paired group states; B: closure controls; C: response separation;
# D: all 100 scaffold errors and median/IQR for each retained context.
# Quantitative grid; layout-only reuse; 183 x 162 mm, Arial, original colors.
# D has a deliberately independent left edge to avoid B's long row-label gutter.
dest <- file.path(OUT,'Fig_6')
svglite::svglite(paste0(dest,'.svg'),width=183/25.4,height=162/25.4,bg='white');print(p);dev.off()
grDevices::cairo_pdf(paste0(dest,'.pdf'),width=183/25.4,height=162/25.4,family='Arial',bg='white');print(p);dev.off()
ragg::agg_png(paste0(dest,'.png'),width=183,height=162,units='mm',res=600,background='white');print(p);dev.off()
ragg::agg_png(paste0(dest,'_preview.png'),width=183,height=162,units='mm',res=300,background='white');print(p);dev.off()
ragg::agg_tiff(paste0(dest,'.tiff'),width=183,height=162,units='mm',res=600,compression='lzw',background='white');print(p);dev.off()
stopifnot(nrow(sep)==100,nrow(err)==300,nrow(ini)==12,nrow(cl)==9)
write_json(list(panel_A_rows=nrow(ini),panel_B_rows=nrow(cl),panel_C_rows=nrow(sep),panel_D_rows=nrow(err),panel_D_summaries=es,source_md5=as.list(tools::md5sum(file.path(D,c('toy_initial_states.csv','closure_summary.csv','same_node_separation.csv','response_errors.csv'))))),paste0(dest,'_data_audit.json'),pretty=TRUE,auto_unbox=TRUE)