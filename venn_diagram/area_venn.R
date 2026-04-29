library(argparser)
library(yaml)

p <- arg_parser("Plot")

p <- add_argument(p, "--inyaml", help = "Input  result file (YAML)", type = "character")

argv <- parse_args(p)

cfg <- yaml::read_yaml(argv$inyaml)

library(ggplot2)
library(ggVennDiagram)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(eulerr)
mat1 = cfg$input$intype1$mat
mat2 = cfg$input$intype2$mat
mat3 = cfg$input$intype3$mat

source('./colormap_utils.R')

print(mat1)
if (! is.null(mat1)){
    df = fread(mat1)
    print(df)
    data = df %>%
  set_names(c("V1", "V2")) %>% 
  group_by(V1) %>%
  summarise(V2 = list(V2), .groups = "drop") %>%
  deframe()

} else if ( ! is.null(mat2)) {
    df = fread(mat2)
    print(df)
    data = df %>% map(~ .x %>%
        discard(~ is.na(.) || . == "") %>%
        unique())
} else if (! is.null(mat3)) {
    df = fread(mat3)
    print(df)
    cutoff = cfg$input$intype3$cutoff
    data <- df %>% rename(gene_col = 1) %>%  
  pivot_longer(-gene_col, names_to = "sample", values_to = "value") %>%
  filter(value > cutoff) %>%
  group_by(sample) %>%
  summarise(gene_id = list(gene_col), .groups = "drop") %>%
  deframe()
}





vd = process_data(Venn(data))
region_df <- venn_region(vd)

col <- get_color_palette(cfg$intersect$fill,length(data))

if (cfg$intersect$border_col == 'follow_fill') {
  edge_col = col
} else {
  edge_col = cfg$intersect$border_col}

if (cfg$intersect$label == 'both') {
  quant_type <- c("percent", "counts")
} else {
  quant_type <- cfg$intersect$label
}
print(edge_col)
print(quant_type)
p <- plot(euler(
     data,
     shape = "circle"),     
     quantities = list(type = quant_type,cex=cfg$intersect$cex,col=cfg$intersect$col),      
     labels=list(cex=cfg$set$cex,col=cfg$set$set_color,fontface=cfg$set$fontface),           
     edges = list(col = edge_col, lwd = cfg$intersect$lwd,lty = cfg$intersect$lty),
     fills = list(fill =col,alpha=cfg$intersect$alpha)
)




ggsave(paste0(cfg$input$out,".png"), plot = p, width = cfg$plot$width, height = cfg$plot$height, dpi = 150, units = "in",bg = "white")


ggsave(paste0(cfg$input$out,".pdf"), plot = p, width = cfg$plot$width, height = cfg$plot$height, units = "in")

region_df$item <- sapply(region_df$item, function(x) paste(x, collapse = ";"))



write.table(region_df,file=paste0(cfg$input$out,".xls"),row.names=F,sep='\t',quote=F)


