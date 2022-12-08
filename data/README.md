# Data description

## Sea level
- **sealevel_full.rds** [Miller et al 2005](https://www.science.org/doi/full/10.1126/science.1116412), taken from the supplement), global mean sea level ranging from 170ma to the recent
- **sealevel_ceno.rds** [Miller et al 2020](https://www.science.org/doi/full/10.1126/sciadv.aaz1346), taken from the supplement), global mean sea level for the Cenozoic


## Temperature 
- **temp_ceno.rds** [Cramer et al 2011](https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2011JC007255&file=jgrc12191-sup-0001-readme.txt), end-Cretaceous deep-ocean temperature
- [Grossman and Joachimski 2022](https://www.nature.com/articles/s41598-022-11493-1#Sec10), supplement 8 with delta18O seawater to calculate surface temperature, whole Phanerozoic
- [Westerhold et al 2020](https://www.science.org/doi/full/10.1126/science.aba6853), raw d180 values of benthic forams, need to be transformed to temperature following Cramer et al 2011, end-Cretaceous to recent deep-ocean temperature
- **temp_full.rds** [Scotese et al 2021](https://www.sciencedirect.com/science/article/pii/S0012825221000027), deep ocean and global average temperatures, whole Phanerozoic, at 2-5 myr resolution
- [Veizher and Prokoph 2015](https://www.sciencedirect.com/science/article/pii/S0012825215000604), d18O of various well-preserved fossils, whole Phanerozoic

## Productivity 
### d13C
- **d13C_ceno.rds** [Westerhold et al 2020](https://www.science.org/doi/full/10.1126/science.aba6853), raw d13C values of benthic forams, sheet S34, end-Cretaceous to recent
- **d13C_full.rds** [Veizher and Prokoph 2015](https://www.sciencedirect.com/science/article/pii/S0012825215000604), raw d13C values of well preserved fossils, whole Phanerozoic
### 87Sr/86Sr
- **SR_full.rds** [McArthur, Howarth, Shields 2012](https://books.google.ch/books?hl=en&lr=&id=1M62_rbq70AC&oi=fnd&pg=PA127&ots=MrPJ841jYQ&sig=gmWJQ8LyKANz3F-rccJAFVt0MiA&redir_esc=y#v=onepage&q&f=false), the lowess fit curve was taken from WebPlotDigitizer, Phanerozoic  
### Diatom richness:
- **diatom_full.rds** downloaded from [Neptune](https://nsb.mfn-berlin.de/) on the 27.11.2022 with TNL taxonomic resolving, without questionable identifications and open-nomenclature taxa and problematic samples or occurrences, using the Gradstein et al al 2012 age scale, filtering out samples with age quality less than medium, and pacman trimming with 5% on both the lower and upper bound

## Outcrop area
- **marine_units_full.rds** Macrostrat (https://macrostrat.org/api/units?age_bottom=150&age_top=0&environ_class=marine&response=longformat=csv) Units with both marine and non-marine environments, need to filter first dat %>% filter(!str_detect(.$environ, "non-marine"))
- **outcrop_full.rds** [Wall, Ivany, Wilkinson 2009](https://www.cambridge.org/core/journals/paleobiology/article/abs/revisiting-raup-exploring-the-influence-of-outcrop-area-on-diversity-in-light-of-modern-samplestandardization-techniques/A1681A7CDCB94EEC7C34161A80E8E4CB) taken from WebPlotDigitizer

## Shelf area
- **cont_area_full.rds** [Kocsis and Scotese 2021](https://www.sciencedirect.com/science/article/pii/S0012825220305092), Phanerozoic flooded continental area as proportion of earths surface area
- **cont_area_ceno.rds** [Miller et al 2005](https://www.science.org/doi/full/10.1126/science.1116412) Cenozoic flooded continental area (10^6 km^2), taken from WebPlotDigitizer
	