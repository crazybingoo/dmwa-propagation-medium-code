suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr);library(readr);library(ggtext);library(svglite);library(ragg);library(jsonlite)})
options(device=function(...) grDevices::cairo_pdf(filename=tempfile(fileext='.pdf'),family='Arial',...))
FIG_TEXT_PT <- 10.5
FIG_PANEL_PT <- 12.5
FIG_HEATMAP_TEXT_PT <- 8
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT/2.845276
INK <- '#202020'; GREY <- '#8A9299'; LIGHT <- '#D3D9DF'
BLUE <- '#0072B2'; ORANGE <- '#D55E00'; GREEN <- '#009E73'; PURPLE <- '#9B73B8'
STAGE_COLORS <- c(Pre='#6E91BF',Early='#83BA7B',Mid='#F1B15E',Late='#D77973',Post='#9E8ABD')
MODEL_COLORS <- c('Full W'=INK,'Coverage-only'=BLUE,'Constant-bias'=ORANGE,'Recovered'=GREEN)
theme_set(theme_classic(base_size=FIG_TEXT_PT,base_family='Arial')+
  theme(axis.text=element_text(size=FIG_TEXT_PT,colour=INK),axis.title=element_text(size=FIG_TEXT_PT,colour=INK),
        axis.title.y=element_text(margin=margin(r=1.5)),axis.title.x=element_text(margin=margin(t=4)),
        axis.line=element_line(linewidth=.35,colour=INK),axis.ticks=element_line(linewidth=.3,colour=INK),
        panel.grid.major=element_line(colour='#EFF2F5',linewidth=.25),panel.grid.minor=element_blank(),
        plot.title=element_markdown(size=FIG_TEXT_PT,face='bold',margin=margin(b=8)),plot.title.position='plot',
        plot.caption=element_text(size=8.4,hjust=0,colour=INK,margin=margin(t=6)),plot.caption.position='plot',
        plot.margin=margin(7,9,7,7),legend.title=element_blank(),legend.text=element_text(size=FIG_TEXT_PT),
        legend.position='bottom',legend.key.width=unit(9,'mm'),legend.key.height=unit(4,'mm'),
        strip.background=element_blank(),strip.text=element_text(size=FIG_TEXT_PT,face='bold')))
ttl <- function(letter,txt) paste0("<span style='font-size:12.5pt'>",letter,"</span>  ",txt)
wrap_figure_notes <- function(p,width=46) {
  if(!is.null(p$labels$caption)) p$labels$caption<-paste(strwrap(p$labels$caption,width=width),collapse='\n')
  if(length(p$patches$plots)) p$patches$plots<-lapply(p$patches$plots,wrap_figure_notes,width=width)
  p
}
