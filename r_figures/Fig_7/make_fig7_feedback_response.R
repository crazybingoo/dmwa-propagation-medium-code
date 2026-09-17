# Fig. 7: existing model-response data, new threshold-led panel order.
# R is the only drawing and export backend. No simulation is run here.
suppressPackageStartupMessages({
 library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
 library(readr); library(ggtext); library(svglite); library(ragg); library(jsonlite)
})
file_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]),winslash="/",mustWork=TRUE))
repo <- normalizePath(file.path(script_dir,"../.."),winslash="/",mustWork=TRUE)
ROOT <- file.path(repo,'mechanism')
D <- Sys.getenv('DMWA_ETA_DATA_DIR',unset=file.path(ROOT,'data/eta'))
OUT <- Sys.getenv('DMWA_FIG7_OUT',unset=file.path(repo,'outputs','Fig_7'))
Q <- OUT;dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
options(device=function(...) grDevices::cairo_pdf(filename=tempfile(fileext='.pdf'),family='Arial',...))
FIG_TEXT_PT <- 10.5; FIG_PANEL_PT <- 12.5; FIG_HEATMAP_TEXT_PT <- 8
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT/2.845276
INK <- '#202020'; GREY <- '#8A9299'; LIGHT <- '#D3D9DF'
BLUE <- '#0072B2'; ORANGE <- '#D55E00'
cols <- c('Chain-attached'=BLUE,'Hub-attached'=ORANGE)
theme_set(theme_classic(base_size=FIG_TEXT_PT,base_family='Arial')+
 theme(axis.text=element_text(size=FIG_TEXT_PT,colour=INK),
 axis.title=element_text(size=FIG_TEXT_PT,colour=INK),
 axis.title.y=element_text(margin=margin(r=1.5)),axis.title.x=element_text(margin=margin(t=4)),
 axis.line=element_line(linewidth=.35,colour=INK),axis.ticks=element_line(linewidth=.3,colour=INK),
 panel.grid.major=element_line(colour='#EFF2F5',linewidth=.25),panel.grid.minor=element_blank(),
 plot.title=element_markdown(size=FIG_TEXT_PT,face='bold',margin=margin(b=8)),plot.title.position='plot',
 plot.subtitle=element_text(size=FIG_TEXT_PT,margin=margin(b=4)),
 plot.margin=margin(7,9,7,7),legend.title=element_blank(),legend.text=element_text(size=FIG_TEXT_PT),
 legend.position='bottom',legend.key.width=unit(9,'mm'),legend.key.height=unit(4,'mm'),
 strip.background=element_blank(),strip.text=element_text(size=FIG_TEXT_PT,face='bold')))
ttl <- function(letter,txt) paste0("<span style='font-size:12.5pt'>",letter,"</span>  ",txt)

par <- read_csv(file.path(D,'example_parameters.csv'),show_col_types=FALSE)
gain_curve <- read_csv(file.path(D,'example_gain_curve.csv'),show_col_types=FALSE)
rates <- read_csv(file.path(D,'example_rates.csv'),show_col_types=FALSE)
curves <- read_csv(file.path(D,'example_trajectories.csv'),show_col_types=FALSE)
pat <- read_csv(file.path(D,'patient_rates.csv'),show_col_types=FALSE)
pm <- read_csv(file.path(D,'patient_matrix_metrics.csv'),show_col_types=FALSE)
stopifnot(nrow(par)==2,nrow(gain_curve)==482,nrow(pat)==360,nrow(pm)==120,
 all(is.finite(gain_curve$direct_eigen_rate)),all(is.finite(curves$log_total_activity)),
 all(is.finite(pat$spectral_rate)),all(is.finite(pat$slope_40_80)),
 max(abs(par$normalized_R-1))<1e-12)
for(mod in unique(gain_curve$model)) stopifnot(all(diff(gain_curve$gain[gain_curve$model==mod])>0))
stopifnot(max(abs(gain_curve$direct_eigen_rate-gain_curve$spectral_rate))<1e-10)

# A: both fixed legal scaffolds; equal mean weight after R normalization.
nodes <- bind_rows(
 data.frame(model='Chain-attached',node=1:5,x=c(.35,1.05,.7,1.62,2.42),y=c(1.14,1.14,.44,.44,.44)),
 data.frame(model='Hub-attached',node=1:5,x=c(4.11,3.55,3.55,4.72,4.72),y=c(.8,1.27,.33,1.27,.33)))
edges <- read_csv(file.path(D,'example_edges.csv'),show_col_types=FALSE)%>%
 left_join(nodes%>%rename(source=node,xs=x,ys=y),by=c('model','source'))%>%
 left_join(nodes%>%rename(target=node,xt=x,yt=y),by=c('model','target'))%>%
 mutate(dx=xt-xs,dy=yt-ys,dist=sqrt(dx^2+dy^2),
 x=xs+.13*dx/dist,y=ys+.13*dy/dist,xend=xt-.13*dx/dist,yend=yt-.13*dy/dist)
labels <- par%>%mutate(x=c(1.27,4.11),y=1.72,
 eta_label=sprintf("italic('η') == '%.3f'",eta))
pA <- ggplot()+
 geom_polygon(data=nodes%>%filter(node<=3),aes(x,y,group=model,fill=model),alpha=.12,colour=NA)+
 geom_segment(data=edges,aes(x,y,xend=xend,yend=yend,colour=model),linewidth=.6)+
 geom_point(data=nodes,aes(x,y,colour=model),size=4.6,shape=21,fill='white',stroke=.6)+
 geom_text(data=nodes,aes(x,y,label=node),size=FIG_HEATMAP_TEXT_PT/2.845276,colour=INK)+
 geom_text(data=labels,aes(x,y,label=model,colour=model),size=FIG_GEOM_TEXT_SIZE,fontface='bold')+
 geom_text(data=labels,aes(x,y=-.03,label=eta_label,colour=model),parse=TRUE,size=FIG_GEOM_TEXT_SIZE)+
 annotate('text',x=2.67,y=-.50,label='5 regions; 5 pairs + 1 triangle each',size=FIG_GEOM_TEXT_SIZE,colour=INK)+
 scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
 coord_cartesian(xlim=c(-.02,5.35),ylim=c(-.68,1.96),clip='off')+
 theme_void(base_family='Arial')+
 theme(legend.position='none',plot.title=element_markdown(size=FIG_TEXT_PT,face='bold',margin=margin(b=8)),
 plot.title.position='plot',plot.subtitle=element_text(size=FIG_TEXT_PT,margin=margin(b=4)),
 plot.margin=margin(7,9,7,7))+
 labs(title=ttl('A','Equal weight, distinct structure'),subtitle='Mean outgoing weight = 1')

# B: all stored gain-curve samples; critical markers are existing analytic roots.
critical <- rates%>%distinct(model,critical_gain)
pB <- ggplot(gain_curve,aes(gain,direct_eigen_rate,colour=model))+
 geom_hline(yintercept=0,colour=GREY,linewidth=.45)+
 geom_segment(data=critical,aes(x=critical_gain,xend=critical_gain,y=-1,yend=0,colour=model),
              inherit.aes=FALSE,linetype='dotted',linewidth=.4)+
 geom_line(linewidth=.7)+
 geom_point(data=critical,aes(x=critical_gain,y=0,colour=model),inherit.aes=FALSE,
            shape=21,fill='white',size=2.5,stroke=.6)+
 annotate('text',x=.18,y=.40,label='Growth',size=FIG_GEOM_TEXT_SIZE,colour=GREY)+
 annotate('text',x=1.06,y=-.86,label='Decay',size=FIG_GEOM_TEXT_SIZE,colour=GREY)+
 scale_colour_manual(values=cols)+
 scale_x_continuous(breaks=c(0,.4,.8,1.2),expand=expansion(mult=c(.02,.025)))+
 scale_y_continuous(breaks=c(-1,-.5,0,.5),expand=expansion(mult=c(.03,.04)))+
 labs(title=ttl('B','Growth and decay thresholds'),
      x=expression('Feedback gain '*italic(g)),y=expression('Dominant growth rate '*italic(s)))+
 theme(legend.position='none')

# C/D retain the original quantitative fields, observations and fitting interval.
curves_plot <- curves%>%mutate(gain=factor(gain,levels=c(.5,.8,.9),labels=c('g = 0.5','g = 0.8','g = 0.9')))
pC <- ggplot(curves_plot,aes(time,log_total_activity,colour=model,linetype=gain))+
 geom_hline(yintercept=0,colour=LIGHT,linewidth=.35)+geom_line(linewidth=.7)+
 scale_colour_manual(values=cols)+scale_linetype_manual(values=c('dotted','dashed','solid'),
 labels=expression(italic(g)==0.5,italic(g)==0.8,italic(g)==0.9))+
 scale_x_continuous(breaks=c(0,20,40,60,80),expand=expansion(mult=c(.01,.02)))+
 labs(title=ttl('C','Response at fixed gains'),x='Time (model units)',y='ln total activity')+
 guides(colour=guide_legend(order=1,nrow=1),linetype=guide_legend(order=2,nrow=1))+
 theme(legend.box='vertical',legend.spacing.y=unit(4,'pt'),legend.margin=margin(0,0,0,0),
       legend.key.width=unit(7,'mm'),legend.text=element_text(size=FIG_TEXT_PT))
lim <- range(c(pat$spectral_rate,pat$slope_40_80))+c(-.025,.025)
pD <- ggplot(pat,aes(spectral_rate,slope_40_80))+
 geom_abline(slope=1,intercept=0,colour=GREY,linewidth=.45)+
 geom_point(shape=21,fill='white',colour=INK,size=1.6,stroke=.4,alpha=.5)+
 scale_x_continuous(limits=lim,breaks=c(-.5,-.3,-.1))+
 scale_y_continuous(limits=lim,breaks=c(-.5,-.3,-.1))+
 labs(title=ttl('D','Model check on patient matrices'),
      subtitle=sprintf('%d matrices; %d model responses',nrow(pm),nrow(pat)),
      x=expression('Analytic growth rate '*italic(s)),y='Numerical late-time growth rate')

# Contract: unchanged model-response evidence; layout-only revision.
# A structural contrast; B thresholds; C responses; D patient-matrix model check.
# Free A from C's axis gutter; collect guides into a balanced full-width footer.
pA <- pA + guides(colour='none',fill='none')
pB <- pB + guides(colour='none')
fig <- wrap_plots(free(pA,side='l'),pB,pC,pD,ncol=2,heights=c(.9,1.1),guides='collect') &
 theme(legend.position='bottom',legend.box='horizontal',legend.spacing.x=unit(12,'pt'))
W_MM <- 183; H_MM <- 165
dest <- file.path(OUT,'Fig_7')
svglite::svglite(paste0(dest,'.svg'),width=W_MM/25.4,height=H_MM/25.4,bg='white');print(fig);dev.off()
grDevices::cairo_pdf(paste0(dest,'.pdf'),width=W_MM/25.4,height=H_MM/25.4,family='Arial',bg='white');print(fig);dev.off()
ragg::agg_png(paste0(dest,'.png'),width=W_MM,height=H_MM,units='mm',res=600,background='white');print(fig);dev.off()
ragg::agg_png(paste0(dest,'_preview.png'),width=W_MM,height=H_MM,units='mm',res=300,background='white');print(fig);dev.off()
ragg::agg_tiff(paste0(dest,'.tiff'),width=W_MM,height=H_MM,units='mm',res=600,compression='lzw',background='white');print(fig);dev.off()
write_json(list(backend='R',width_mm=W_MM,height_mm=H_MM,dpi=600,
 source_counts=list(structures=nrow(par),gain_rows=nrow(gain_curve),trajectory_rows=nrow(curves),
 matrices=nrow(pm),patient_responses=nrow(pat)),excluded_rows=0,
 maximum_gain_curve_identity_error=max(abs(gain_curve$direct_eigen_rate-gain_curve$spectral_rate)),
 maximum_patient_slope_error=max(pat$slope_abs_error),files=paste0(dest,c('.svg','.pdf','.png','.tiff'))),
 file.path(Q,'Fig_7_data_export_audit.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
