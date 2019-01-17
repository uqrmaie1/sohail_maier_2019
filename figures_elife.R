#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Code for reproducing the figures in
# "Polygenic adaptation on height is overestimated due to uncorrected stratification in genome-wide association studies"
#
# Robert Maier
#
# 2019/01/16
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(magrittr)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 1b
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat2sds_f1 = read_delim('data/figure1_d1.tsv', '\t')
allbj2 = read_delim('data/figure1_d2.tsv', '\t')

fc = function(x) formatC(x, format = 'e', digits = 0)

p1_2 = ggplot(dat2sds_f1, aes(spbin, mntSDS)) +
  xlab('P value bin (high to low P values)') +
  ylab('Height increasing tSDS\n(bin average)') +
  geom_hline(yintercept=0, col='grey') +
  geom_point(size=.25) +
  facet_grid(. ~ nam4, scales='free_x', labeller=as_labeller(function(s) gsub('.+', '', s))) +
  geom_smooth(method='lm', se=F, aes(col=nam4)) +
  theme(legend.position='none', panel.background = element_blank(), panel.grid = element_blank(),
        panel.spacing = unit(0, 'lines'), axis.line = element_line(), strip.background=element_blank()) +
  scale_x_reverse(breaks=0:1) +
  scale_alpha_manual(values=c(.5, 1)) +
  geom_text(aes(x=1, y=.5, label=paste0('Spearman r: ', round(average, 3), '\nP-value: ', fc(pvalue))), data=allbj2, hjust=0)

pdf('figure1b.pdf', width=5, height=4)
print(p1_2)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Figure 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat2pc = read_table2('data/figure2_d1.tsv')
print(load('data/figure2_d2.RData'))

numpcs = 20
bluered = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", 
                                 "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                 "#4393C3", "#2166AC", "#053061")))


dat2pc_beta = filter(dat2pc, nam1 != 'SDS')

p1 = ggplot(dat2pc_beta, aes(PC, corr, fill=gbrtsicorr)) +
  ylab(expression(rho(PC, beta))) +
  scale_fill_gradientn(expression(rho(PC, Delta[AF])), colors=bluered(20),
                       limits=c(-max(abs(dat2pc_beta$gbrtsicorr), na.rm=T),
                                max(abs(dat2pc_beta$gbrtsicorr), na.rm=T)))

dat2pc_sds = filter(dat2pc, nam1 == 'SDS')

p2 = ggplot(dat2pc_sds, aes(PC, corr, fill=gbrtsicorr)) +
  ylab(expression(rho(PC, SDS))) +
  scale_fill_gradientn(expression(rho(PC, Delta[AF])), colors=bluered(20),
                       limits=c(-max(abs(dat2pc_sds$gbrtsicorr), na.rm=T), max(abs(dat2pc_sds$gbrtsicorr), na.rm=T)))

comm1 = list(geom_bar(stat='identity', position='dodge', col='black'),
             facet_grid(. ~ nam2, scales='free_y', labeller = label_parsed),
             geom_errorbar(aes(ymin=corr-1.96*corse, ymax=corr+1.96*corse), width=.1),
             geom_text(aes(x=PC, y=max(corr), label=ifelse(pval < .05, ifelse(pval < 0.05/numpcs, '*', '°'), ''))),
             xlab('PC'), scale_x_continuous(breaks=seq(1,numpcs,2), expand=c(0,0.1)),
             theme(panel.grid = element_blank(), panel.background = element_blank(),
                   axis.line.x = element_line(), strip.background = element_blank()),
             geom_vline(xintercept=0))

# heatmap
bins = 30

p3 = ggplot(filter(datcomb, nam1 != 'SDS', !is.na(gbrbin))) +
  scale_fill_gradient2(name=expression(bar(beta)), na.value='white',
                       high = scales::muted("red"), low = scales::muted("blue"), mid='#FAFAFA')

p4 = ggplot(filter(datcomb, nam1 == 'SDS', !is.na(gbrbin))) +
  scale_fill_gradient2(name=expression(bar(SDS)), na.value='white',
                       high = scales::muted("red"), low = scales::muted("blue"), mid='#FAFAFA')

brks = sort(unique(datcomb$tsi_binmean))[c(1,seq(3,bins,5))]

comm2 = list(geom_tile(aes(tsi_binmean, gbr_binmean, fill=mb_scaled_na)),
             facet_grid(. ~ nam2, labeller = label_parsed),
             theme(axis.text.x = element_text(angle = 90), panel.grid = element_blank(),
                   axis.text=element_text(size=8), panel.background = element_blank(),
                   axis.line.x = element_line(), strip.background = element_blank()),
             ylab('MAF bin GBR'), xlab('MAF bin TSI'),
             geom_abline(linetype='dotted', col='black'),
             geom_vline(xintercept=c(0.5, 1.5)),
             geom_hline(yintercept=1.5),
             scale_x_discrete(breaks=brks),
             scale_y_discrete(breaks=brks))

w1 = 0.62; w2 = 1-w1
g1 = grid.arrange(ggarrange(ggplotGrob(p1 + comm1), labels='a'), widths=w1)
g2 = grid.arrange(ggarrange(ggplotGrob(p2 + comm1), labels='b'), widths=w2)
g3 = grid.arrange(ggarrange(ggplotGrob(p3 + comm2), labels='c'), widths=w1)
g4 = grid.arrange(ggarrange(ggplotGrob(p4 + comm2), labels='d'), widths=w2)
gr1 = cbind(g1, g2, size = "first")
gr2 = cbind(g3, g4, size = "first")
g = rbind(gr1, gr2, size = "first")
g$widths = unit.pmax(gr1$widths, gr2$widths)

pdf('figure2.pdf', width=11, height=7)
grid.newpage()
grid.draw(g)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat2sds = read_delim('data/figure3_d1.tsv', '\t')
dat3sds = read_delim('data/figure3_d2.tsv', '\t')
sums_comb_fig3b = read_delim('data/figure3_d3.tsv', '\t')

# Figure 3a

p1a = ggplot(dat2sds, aes(spbin, mntSDS)) +
  geom_text(aes(x=1, y=max(dat2sds$mntSDS),
                label=paste0('beta:~', fc(slope_tSDS)), alpha=pslope_tSDS < .05),
            data=dat3sds, hjust=0, vjust=1, parse=T, size=3) +
  geom_text(aes(x=1, y=max(dat2sds$mntSDS)-.1*diff(range(dat2sds$mntSDS)),
                label=paste0('P:~',fc(pslope_tSDS)), alpha=pslope_tSDS < .05),
            data=dat3sds, hjust=0, vjust=1, parse=T, size=3) +
  xlab('') +
  ylab(expression(atop(tSDS,of~height~increasing~allele)))

p2a = ggplot(dat2sds, aes(spbin, mndiffsign)) +
  geom_text(aes(x=1, y=max(dat2sds$mndiffsign),
                label=paste0('beta:~', fc(slope_afdiff)), alpha=pslope_afdiff < .05),
            data=dat3sds, hjust=0, vjust=1, parse=T, size=3) +
  geom_text(aes(x=1, y=max(dat2sds$mndiffsign)-.1*diff(range(dat2sds$mndiffsign)),
                label=paste0('P:~',fc(pslope_afdiff)), alpha=pslope_afdiff < .05),
            data=dat3sds, hjust=0, vjust=1, parse=T, size=3) +
  xlab('P value bin (high to low P-values)') +
  ylab(expression(atop(Delta[AF]~'(GBR - TSI)',of~height~increasing~allele)))

comma = list(geom_hline(yintercept=0, col='grey'),
             geom_point(size=.25),
             facet_grid(. ~ nam4, scales='free_x',
                        labeller=as_labeller(function(s) gsub('[0-9]+\\. +', '', s))),
             geom_smooth(method='lm', se=F, aes(col=nam4)),
             theme(legend.position='none', panel.background = element_blank(),
                   panel.grid = element_blank(), panel.spacing = unit(0, 'lines'),
                   axis.line = element_line(), strip.background=element_blank()),
             scale_x_reverse(breaks=0:1), scale_alpha_manual(values=c(1, 1)))


g1 = grid.arrange(ggarrange(ggplotGrob(p1a + comma), labels='a'))
g2 = grid.arrange(ggarrange(ggplotGrob(p2a + comma), labels='c'))
g = rbind(g1, g2, size = "first")
g$widths = unit.pmax(g1$widths, g2$widths)


# Figure 3b

lwd = .5
col1 = 'darkgoldenrod2'
col2 = 'grey'

p1b = sums_comb_fig3b %>% ggplot(aes(tSDS)) +
  stat_function(fun=dnorm, lwd=lwd, col=col2) +
  geom_vline(xintercept=c(0, mean(sums_comb_fig3b$tSDS, na.rm=T)), col=c(col2, col1), lwd=lwd, linetype='dashed') +
  scale_x_continuous(expression(tSDS~'('~UKB~sibs~WB~beta~')'), limits=c(-4,4)) +
  annotate('text', size=2.5, x=1, y=.4, hjust=0, vjust=1,
           label=paste0('n=', nrow(sums_comb_fig3b),'\nmean=',round(mean(sums_comb_fig3b$tSDS), 2),
                        '\np=',fc(t.test(sums_comb_fig3b$tSDS)$p.value)))

p2b = sums_comb_fig3b %>% ggplot(aes(afdiff)) +
  geom_vline(xintercept=c(0, mean(sums_comb_fig3b$afdiff, na.rm=T)), col=c(col2, col1), lwd=lwd, linetype='dashed') +
  scale_x_continuous(expression(Delta[AF]~'(GBR - TSI)')) +
  annotate('text', size=2.5, x=.1, y=5, hjust=0, vjust=1,
           label=paste0('n=', sum(!is.na(sums_comb_fig3b$afdiff)),'\nmean=',
                        round(mean(sums_comb_fig3b$afdiff, na.rm=T), 3),
                        '\np=',formatC(wilcox.test(sums_comb_fig3b$afdiff)$p.value, format = 'g', digits = 2)))

commb = list(geom_density(lwd=lwd, col=col1),
             ylab('density'),
             theme(panel.grid = element_blank(), plot.title=element_text(size=10),
                   axis.title.x = element_text(size=10), panel.background = element_blank(),
                   axis.line = element_line()),
             ggtitle(paste0('UKB p < 5e-8, LD-pruned')))

pdf('figure3.pdf', width=9, height=5)
grid.arrange(arrangeGrob(g,
                         ggarrange(ggplotGrob(p1b + commb), labels='b'),
                         ggarrange(ggplotGrob(p2b + commb), labels='d'),
                         layout_matrix = matrix(c(1,1,2,3),2), widths=c(.75, .25)))
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Figure 4: POPRES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat_fig4b = read_delim('data/figure4_d1.tsv', '\t')
dat3 = read_delim('data/figure4_d2.tsv', '\t')
dat4 = read_delim('data/figure4_d3.tsv', '\t')

gg_color_hue = function(n, l=65, c=100) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=c)[1:n]
}

set4 = c('giant', 'ukb_neale1', 'ukbfam_wb', 'ukbsparse_all_nopcs')

f4b = dat_fig4b %>%
  filter(nam1 %in% set4, clump %in% c('UKB sig', 'clump 0.01')) %>%
  mutate(clump=recode(clump, 'UKB sig'='p<5%*%10^{-8}~(UKB)', 'clump 0.01'='p<0.01'),
         country2=recode(country, 'Italy'=' Italy', 'United Kingdom'=' United Kingdom', 'Spain'=' Spain', .default='other')) %>%
  ggplot(aes(latitude, mn)) +
  geom_smooth(method = 'lm', col='grey', se=F) +
  geom_point(aes(color=country2)) +
  geom_errorbar(width=.01, aes(ymin=mn-1.96*se2, ymax=mn+1.96*se2, col=country2)) +
  facet_grid(clump ~ nam_set4, scales='fixed', labeller=label_parsed) +
  geom_text(parse=T, hjust=0, vjust=1,
            aes(x=39, y=max(top)*.95, label=paste0('P[Qx]==', ifelse(qx<1e-250, '0', paste0(gsub('e', '%*%10^{', fc(qx)), '}')))),
            data=dat3) +
  geom_text(parse=T, hjust=0, vjust=1,
            aes(x=39, y=max(top)*.7, label=paste0('P[lat]==', gsub('e', '%*%10^{', fc(lat)), '}')),
            data=dat3) +
  xlab(paste0('Latitude (',format(min(dat4$latitude), nsmall=1), ' to ', format(max(dat4$latitude), nsmall=1), ')' )) +
  ylab('Polygenic height score') +
  theme_bw() +
  scale_x_continuous(breaks=dat4$latrank2, labels=dat4$country) +
  scale_y_continuous(expand=c(0,0)) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1, size=5),
        legend.position = 'bottom', panel.spacing = unit(0, 'lines'),
        panel.border = element_blank(), axis.line.y = element_line(), strip.background = element_blank()) +
  scale_color_manual('', values=c(gg_color_hue(3)[1], gg_color_hue(3)[2], gg_color_hue(3)[3], 'black')) +
  geom_hline(yintercept = min(dat3$bottom))

ggsave('figure4.pdf', f4b, height=5, width=8.5)

