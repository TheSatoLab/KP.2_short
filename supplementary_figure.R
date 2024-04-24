library(tidyverse)
library(stats)
library(base)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(datasets)
library(grDevices)
library(patchwork)
library(RColorBrewer)
library(cmdstanr)


set_cmdstan_path("cmdstan")
cmdstan_path()

########## args ##########
#Change when using new input
download_date <- "2024-04-11"
out_prefix <- "2024_04_11"
date_w_space <- "20240411"

##input
metadata.name <- paste(out_prefix,"/metadata_tsv_",out_prefix,"/metadata.tsv",sep = "")
mut.info.name <- paste(out_prefix,"/metadata_tsv_",out_prefix,"/metadata.mut_long.tsv",sep = "")
stan_f.name <- ""
#output
metadata_500.name <- paste(out_prefix,"/metadata_500.tsv",sep = "")
mutation_figure.name <- paste(out_prefix,'.KP.2.pdf', sep = "")

dir <- paste("output/",out_prefix,sep = "")
setwd(dir)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########parameters##########
##general
core.num <- 4
variant.ref <- "JN.1"

##period to be analyzed
date.start <- as.Date("2023-08-01")
date.end <- as.Date(download_date)

##min numbers
limit.count.analyzed <- 50

##Transmissibility
bin.size <- 1
generation_time <- 2.1

##model
multi_nomial_model <- cmdstan_model(stan_f.name)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########## Filtering metadata #########
## metadata
metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)

# i) only ‘original passage’ sequences
# iii) host labelled as ‘Human’
# iv) sequence length above 28,000 base pairs
# v) proportion of ambiguous bases below 2%.

metadata.filtered <- metadata %>%
  distinct(Accession.ID,.keep_all=T) %>%
  filter(Host == "Human",
         !N.Content > 0.02 | is.na(N.Content),
         str_length(Collection.date) == 10,
         Sequence.length > 28000,
         Passage.details.history == "Original",
         Pango.lineage != "",
         Pango.lineage != "None",
         Pango.lineage != "Unassigned",
         !str_detect(Additional.location.information,"[Qq]uarantine")
  )

metadata.filtered <- metadata.filtered %>%
  mutate(Collection.date = as.Date(Collection.date),
         region = str_split(Location," / ",simplify = T)[,1],
         country = str_split(Location," / ",simplify = T)[,2],
         state = str_split(Location," / ",simplify = T)[,3])

metadata.filtered <- metadata.filtered[!duplicated(metadata.filtered$Virus.name),]

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########## Mutation frequency plot for tree reconstruction #########
## mut.info
## make sure to use python script that includes all needed period to make mut.info
metadata.filtered.2 <- metadata.filtered %>% filter(country == "USA")

lineage.interest.v <- c("BA.2.86.1", "JN.1","KP.2")

mut.info <- fread(mut.info.name,header=T,sep="\t",quote="",check.names=T)

mut.info.merged <- mut.info %>% inner_join(metadata.filtered.2 %>% select(Id=Accession.ID,Pango.lineage),by="Id") %>% mutate(mut = str_replace(mut, ":", "_")) %>% mutate(prot=str_split(mut, "_", simplify=T)[,1], mut.mod=gsub("[A-Z]", "", str_split(mut, "_", simplify=T)[,2], ignore.case=TRUE))
mut.info.merged <- mut.info.merged %>% as.data.frame() %>% filter(Pango.lineage %in% lineage.interest.v) #%>% slice_sample(n = 5000)

metadata.filtered.2 <- metadata.filtered.2 %>% filter(Accession.ID %in% as.character(mut.info.merged$Id))
metadata.filtered.2 <- metadata.filtered.2 %>% filter(!(Pango.lineage=="XBB.1" & ((! grepl("Spike_F486S",AA.Substitutions)) | (!grepl("Spike_T478K",AA.Substitutions)))))
metadata.filtered.2 <- metadata.filtered.2 %>% distinct(Accession.ID,.keep_all=T)

liste_mut.name <- paste(date_w_space,"_list_mut.tsv",sep = "")
metadata.filtered.epi_mut_set.list <- paste(metadata.filtered.2$Accession.ID)
write.table(metadata.filtered.epi_mut_set.list, liste_mut.name, col.names=F, row.names=F, sep="\n", quote=F)

mut.info.merged <- mut.info.merged %>% filter(Id %in% as.character(metadata.filtered.2$Accession.ID))

count.pango.df <- metadata.filtered.2 %>% group_by(Pango.lineage) %>% summarize(count.pango = n())
count.pango_mut.df <- mut.info.merged %>% group_by(Pango.lineage,mut) %>% summarize(count.pango_mut = n())

count.pango_mut.df.merged <- count.pango_mut.df %>% inner_join(count.pango.df,by="Pango.lineage")
count.pango_mut.df.merged <- count.pango_mut.df.merged %>% mutate(mut.freq = count.pango_mut / count.pango)

mut.interest.v <- count.pango_mut.df.merged %>% filter(mut.freq > 0.5) %>% pull(mut) %>% unique()

count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged %>% filter(mut %in% mut.interest.v)
count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% mutate(mut.freq.binary = ifelse(mut.freq > 0.2,1,0))

count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% mutate(gene = gsub("_[^_]+","",mut), AA_change = gsub("[^_]+_","",mut), pos = gsub("[A-Za-z]","",AA_change) %>% as.numeric())

mut.spread <- count.pango_mut.df.merged.filtered %>% select(Pango.lineage,mut,gene,AA_change,mut.freq.binary) %>% spread(key = Pango.lineage, value = mut.freq.binary)
mut.spread[is.na(mut.spread)] <- 0

mut.interest.v2 <- mut.spread[apply(mut.spread[,4:ncol(mut.spread)],1,sd) > 0,]$mut

count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% filter(mut %in% mut.interest.v2)

count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>%
  mutate(mut_mod = ifelse(pos %in% 69:70,"HV69-70del",as.character(AA_change)))


count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% group_by(Pango.lineage,gene,mut_mod) %>% summarize(mut.freq = mean(mut.freq)) %>% ungroup() %>% mutate(mut = paste(gene,mut_mod,sep="_"))

mut.order.v <- count.pango_mut.df.merged.filtered %>% select(gene,mut_mod,mut) %>% mutate(pos = gsub("\\-.+","",gsub("[A-Za-z]","",mut_mod))) %>% distinct(mut_mod,.keep_all =T) %>% arrange(gene,as.numeric(pos)) %>% pull(mut)

count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% mutate(mut = factor(mut,levels=rev(mut.order.v)))


count.pango_mut.df.merged.filtered <- count.pango_mut.df.merged.filtered %>% mutate(Pango.lineage = factor(Pango.lineage,levels=lineage.interest.v))

g <- ggplot(count.pango_mut.df.merged.filtered, aes(x = Pango.lineage, y = mut, fill = mut.freq))
g <- g + geom_tile()
g <- g + scale_fill_gradientn(colours=brewer.pal(9, "BuPu"),limits=c(0,1))
g <- g + theme_set(theme_classic(base_size = 10, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()
)
g <- g + theme(
  legend.key.size = unit(0.3, 'cm'), #change legend key size
  legend.key.height = unit(0.3, 'cm'), #change legend key height
  legend.key.width = unit(0.3, 'cm'), #change legend key width
  legend.title = element_text(size=6), #change legend title font size
  legend.text = element_text(size=6))
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + labs(fill = "Freq.")
g <- g + xlab("") + ylab("")
g <- g + theme(axis.ticks=element_line(colour = "black"),
               axis.text=element_text(colour = "black"))
g

pdf(mutation_figure.name,width=2.5,height=2.5)
plot(g)
dev.off()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###mutate a new Pango.lineage column that sets "JN.1+R346T" and "JN.1+F456L"
metadata.filtered.2 <- metadata.filtered %>% 
  mutate(
    new_Pango_lineage = ifelse(
      Variant == "VOI GRA (JN.1+JN.1.*) first detected in Luxembourg/Iceland" & 
        grepl("Spike_R346T", AA.Substitutions) & 
        !grepl("Spike_F456L", AA.Substitutions),
      "JN.1+R346T",
      ifelse(
        Variant == "VOI GRA (JN.1+JN.1.*) first detected in Luxembourg/Iceland" & 
          grepl("Spike_F456L", AA.Substitutions) & 
          !grepl("Spike_R346T", AA.Substitutions),
        "JN.1+F456L",
      ifelse(
        Variant == "VOI GRA (JN.1+JN.1.*) first detected in Luxembourg/Iceland" & 
          !grepl("Spike_F456L", AA.Substitutions) & 
          !grepl("Spike_R346T", AA.Substitutions),
        "JN.1",
        Pango.lineage
      ))
    )
  )

metadata.filtered.2 <- metadata.filtered.2 %>% filter(Collection.date >= date.start, Collection.date <= date.end)
head(metadata.filtered.2)

filter(metadata.filtered.2, new_Pango_lineage=="JN.1+R346T")
JN1_R346T_variants <- metadata.filtered.2[metadata.filtered.2$new_Pango_lineage == "JN.1+R346T", "Pango.lineage"]
unique(JN1_R346T_variants)

filter(metadata.filtered.2, new_Pango_lineage=="JN.1+F456L")
JN1_F456L_variants <- metadata.filtered.2[metadata.filtered.2$new_Pango_lineage == "JN.1+F456L", "Pango.lineage"]
unique(JN1_F456L_variants)


########## Predict effective reproductive number (Re) ##########
## modeling
lineage.interest.v <- c("HK.3","BA.2.86.1", "JN.1", "JN.1+R346T","JN.1+F456L","KP.2")
country.interest <- c("USA", "United Kingdom", "Canada")
plot.l <- list()

for (a in country.interest) {

metadata.filtered.interest <- metadata.filtered.2 %>% filter(country == a)
metadata.filtered.interest <- metadata.filtered.interest %>% mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1, date.bin = cut(date.num,seq(0,max(date.num),bin.size)), date.bin.num = as.numeric(date.bin))
metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin))

##filter variants more than the set limit.count.analyzed
count.interest <- metadata.filtered.interest %>% group_by(new_Pango_lineage) %>% summarize(count = n())
count.interest.50 <- count.interest %>% filter(count >= limit.count.analyzed | new_Pango_lineage=="KP.2" | new_Pango_lineage=="JN.1+F456L" | new_Pango_lineage=="JN.1+R346T")
variant.more.50.v <- unique(c(count.interest.50$new_Pango_lineage, variant.ref)) 
metadata.filtered.interest <- metadata.filtered.interest %>% filter(new_Pango_lineage %in% variant.more.50.v)

list_Re.name <- paste(date_w_space, "_list_Re_",a,".tsv",sep = "")
metadata.filtered.epi_set.list <- paste(metadata.filtered.interest$Accession.ID)
write.table(metadata.filtered.epi_set.list, list_Re.name, col.names=F, row.names=F, sep="\n", quote=F)

##count variants per day
metadata.filtered.interest.bin <- metadata.filtered.interest %>% group_by(date.bin.num, new_Pango_lineage) %>% summarize(count = n()) %>% ungroup()

metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin %>% spread(key=new_Pango_lineage,value = count)
metadata.filtered.interest.bin.spread[is.na(metadata.filtered.interest.bin.spread)] <- 0
metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin.spread

X <- as.matrix(data.frame(X0 = 1, X1 = metadata.filtered.interest.bin.spread$date.bin.num))

Y <- metadata.filtered.interest.bin.spread %>% select(- date.bin.num)
Y <- Y[,c(variant.ref,colnames(Y)[-which(colnames(Y)==variant.ref)])]

count.group <- apply(Y,2,sum)
count.total <- sum(count.group)
prop.group <- count.group / count.total

Y <- Y %>% as.matrix()
apply(Y,2,sum)

group.df <- data.frame(group_Id = 1:ncol(Y), group = colnames(Y))

Y_sum.v <- apply(Y,1,sum)

data.stan <- list(K = ncol(Y),
                  D = 2,
                  N = nrow(Y),
                  X = X,
                  Y = Y,
                  generation_time = generation_time,
                  bin_size = bin.size,
                  Y_sum = c(Y_sum.v))

fit.stan <- multi_nomial_model$sample(
  data=data.stan,
  iter_sampling=4000,
  iter_warmup=1000,
  seed=1234,
  parallel_chains = 4,
  #adapt_delta = 0.99,
  max_treedepth = 15,
  #pars=c('b_raw'),
  chains=4)

#growth rate
stat.info <- fit.stan$summary("growth_rate") %>% as.data.frame()
stat.info$Nextclade_pango <- colnames(Y)[2:ncol(Y)]

stat.info.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame() %>% rename(q2.5 = `2.5%`, q97.5 = `97.5%`)
stat.info <- stat.info %>% inner_join(stat.info.q,by="variable")

out.name <- paste('growth_rate.2023-04-10.wo_strata.', a, '.txt', sep='')
write.table(stat.info, out.name, col.names=T, row.names=F, sep="\t", quote=F)


draw.df.growth_rate <- fit.stan$draws("growth_rate", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.growth_rate.long <- draw.df.growth_rate %>% gather(key = class, value = value)

draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group_Id = str_match(draw.df.growth_rate.long$class,'growth_rate\\[([0-9]+)\\]')[,2] %>% as.numeric() + 1)
draw.df.growth_rate.long <- merge(draw.df.growth_rate.long,group.df,by="group_Id") %>% select(value,group)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% group_by(group) %>% filter(value>=quantile(value,0.005),value<=quantile(value,0.995))
draw.df.growth_rate.long <- rbind(data.frame(group=variant.ref,value=1),draw.df.growth_rate.long)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% filter(group %in% lineage.interest.v)

draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group = factor(group, levels=lineage.interest.v))

col.v <- brewer.pal(length(lineage.interest.v) + 1, "Set1")[c(1:6)]
g1 <- ggplot(draw.df.growth_rate.long,aes(x=group,y=value,color=group,fill=group))
g1 <- g1 + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
g1 <- g1 + geom_violin(alpha=0.6,scale="width")
g1 <- g1 + stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975), size=0.5,fatten =1.5)
# g1 <- g1 + scale_fill_manual(valus = color_maping)
g1 <- g1 + scale_color_manual(values=col.v)
g1 <- g1 + scale_fill_manual(values=col.v)
g1 <- g1 + theme_classic()
g1 <- g1 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8))
g1 <- g1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g1 <- g1 + ggtitle(a)
g1 <- g1 + xlab('') + ylab('Relative Re')
g1 <- g1 + theme(legend.position = 'none')
g1 <- g1 + scale_y_continuous(limits=c(0.8,1.8),breaks=c(0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8))
g1

pdf.name <- paste('Re_violin_plot.wo_strata.', a, '.pdf', sep='')
plot.l[[pdf.name]] <- g1
pdf(pdf.name,width=3.5,height=5.0)
plot(g1)
dev.off()

data_Id.df <- data.frame(data_Id = 1:length(X), date_Id = X[,2], Y_sum = Y_sum.v, date = as.Date(X[,2],origin=date.start)-1)


data.freq <- metadata.filtered.interest.bin %>% rename(group = new_Pango_lineage) %>% group_by(date.bin.num) %>% mutate(freq = count / sum(count))
data.freq <- data.freq %>% mutate(date = as.Date(date.bin.num,origin=date.start)-1)


draw.df.theta <- fit.stan$draws("theta", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.theta.long <- draw.df.theta %>% gather(key = class, value = value)
draw.df.theta.long <- draw.df.theta.long %>% mutate(data_Id = str_match(class,'theta\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(),
                                                    group_Id = str_match(class,'theta\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric())

draw.df.theta.long <- draw.df.theta.long %>% inner_join(data_Id.df %>% select(data_Id,date), by = "data_Id")

draw.df.theta.long.sum <- draw.df.theta.long %>% group_by(group_Id, date) %>% summarize(mean = mean(value),ymin = quantile(value,0.025),ymax = quantile(value,0.975))
draw.df.theta.long.sum <- draw.df.theta.long.sum %>% inner_join(group.df,by="group_Id")

draw.df.theta.long.sum.filtered <- draw.df.theta.long.sum %>% filter(group %in% lineage.interest.v) %>% mutate(group = factor(group,levels=lineage.interest.v))

g2 <- ggplot(draw.df.theta.long.sum.filtered,aes(x=date, y = mean, fill=group, color = group))
g2 <- g2 + geom_ribbon(aes(ymin=ymin,ymax=ymax), color=NA,alpha=0.2)
g2 <- g2 + geom_line(linewidth=0.3)
g2 <- g2 + scale_x_date(date_labels = "%y-%m", date_breaks = "1 months", date_minor_breaks = "1 month", limits = c(date.start, date.end))
g2 <- g2 + theme_classic()
g2 <- g2 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8)
)
g2 <- g2 + ggtitle(a)
g2 <- g2 + scale_color_manual(values = col.v)
g2 <- g2 + scale_fill_manual(values = col.v)
g2 <- g2 + scale_y_continuous(limits=c(0,1.0),breaks=c(0,0.25,0.5,0.75,1.0))
g2

pdf.name <- paste('lineage_dynamics.wo_strata.', a, '.pdf', sep='')
plot.l[[pdf.name]] <- g2
pdf(pdf.name,width=5,height=5)
plot(g2)
dev.off()
}

# pdf("Re_violin_plot.wo_strata.all.pdf", width=30, height=20)
rep_num_plot <- plot.l[[1]] + plot.l[[3]] + plot.l[[5]] + plot.l[[7]] + plot_layout(ncol=1)
# rep_num_plot
# dev.off()

# pdf("lineage_dynamics.wo_strata.all.pdf", width=30,height=20)
epi_dynamics_plot <- plot.l[[2]] + plot.l[[4]] + plot.l[[6]] + plot.l[[8]] + plot_layout(nrow=2, ncol=2)
# epi_dynamics_plot 
# dev.off()

# wrap_plots(plot.l,ncol=1)
