#################### setting colors for project #############



# updating after new polariation for ms and mm (separating) and pb and bb (with pb mapped to pb)
niceLabels= data.frame(label=c("brown_bear_ABC","brown_bear_EUR","polar_bear_PB","fin_whale_ENP","fin_whale_GOC","Gorilla_gorilla","humans_AFR","humans_AMR","humans_EAS","humans_EUR","humans_SAS","Mus_musculus_Mmc","Mus_musculus_Mmd","Mus_musculus_Mmm","Mus_spretus_Ms","Pan_paniscus","Pan_troglodytes", "Pongo_abelii","Pongo_pygmaeus","vaquita","wolves","ucla_wolves"),niceLabel=c("brown bear (ABC)","brown bear (EUR)","polar bear","fin whale (ENP)","fin whale (GOC)","gorilla","human (AFR)","human (AMR)","human (EAS)","human (EUR)","human (SAS)","mouse (Mmc)","mouse (Mmd)","mouse (Mmm)", "mouse (Ms)","bonobo","chimpanzee","orangutan (Sumatran)","orangutan (Bornean)","vaquita","wolf (Broad)","wolf (UCLA)"),shortSpeciesLabel=c("brown bear","brown bear","polar bear","fin whale","fin whale","gorilla","human","human","human","human","human","house mouse","house mouse","house mouse","Algerian mouse","bonobo","chimpanzee","S. orangutan","B. orangutan","vaquita",'wolf (Broad)',"wolf (UCLA)"))

#colorList=list("brown bear (ABC)"="chocolate4","brown bear (EUR)"="chocolate4","polar bear"="thistle","fin whale (ENP)"="slateblue","fin whale (GOC)"="slateblue","gorilla"="magenta3","human (AFR)"="darkcyan","human (AMR)"="darkcyan","human (EAS)"="darkcyan","human (EUR)"="darkcyan","human (SAS)"="darkcyan","mouse (Mmc)"="orangered","mouse (Mmd)"="orangered","mouse (Mmm)"="orangered", "mouse (Ms)"="orange2","bonobo"="gold3","chimpanzee"="lightsalmon","orangutan (Sumatran)"="darkgoldenrod2","orangutan (Bornean)"="olivedrab3","vaquita"="cyan3")


#colorList=list("brown bear (ABC)"="darkorchid3","brown bear (EUR)"="darkorchid3","polar bear"="thistle","fin whale (ENP)"="slateblue","fin whale (GOC)"="slateblue","gorilla"="magenta3","human (AFR)"="tomato2","human (AMR)"="tomato2","human (EAS)"="tomato2","human (EUR)"="tomato2","human (SAS)"="tomato2","mouse (Mmc)"="bisque4","mouse (Mmd)"="bisque4","mouse (Mmm)"="bisque4", "mouse (Ms)"="tan3","bonobo"="gold3","chimpanzee"="lightcoral","orangutan (Sumatran)"="deeppink1","orangutan (Bornean)"="palevioletred2","vaquita"="cyan3")


# idea: make apes sunset colors; whales; bears purples; mice browns?

# make a pie chart:
#uniquecolors=unique(unlist(unname(colorList)))
#area <- rep(1,length(uniquecolors))
#pie(area,col=uniquecolors)

#display.brewer.all(colorblindFriendly = TRUE) # thsee are color blind friendly(?)
#n <- 12
#col=brewer.pal(n=n,"Paired")
#area <- rep(1,n)
#pie(area, col = col)

#olorList=list("brown bear (ABC)"=col[12],"brown bear (EUR)"=col[12],"polar bear"="thistle","fin whale (ENP)"=col[2],"fin whale (GOC)"=col[2],"gorilla"=col[11],"human (AFR)"=col[9],"human (AMR)"=col[9],"human (EAS)"=col[9],"human (EUR)"=col[9],"human (SAS)"=col[9],"mouse (Mmc)"=col[4],"mouse (Mmd)"=col[4],"mouse (Mmm)"=col[4], "mouse (Ms)"=col[3],"bonobo"=col[5],"chimpanzee"=col[6],"orangutan (Sumatran)"=col[7],"orangutan (Bornean)"=col[8],"vaquita"=col[1])
 # adding in thistle

#http://mkweb.bcgsc.ca/biovis2012/
#http://mkweb.bcgsc.ca/colorblind/palettes.mhtml#page-container
### 15 color palette : 
col <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")

pie(rep(1,15), col=col)
# yellow is too bright so sub it out :
col <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
           "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
           "#920000","#924900","#db6d00","#24ff24","yellow3")
colorList=list("brown bear (ABC)"=col[12],"brown bear (EUR)"=col[12],"polar bear"=col[2],"fin whale (ENP)"=col[7],"fin whale (GOC)"=col[7],"gorilla"=col[15],"human (AFR)"=col[14],"human (AMR)"=col[14],"human (EAS)"=col[14],"human (EUR)"=col[14],"human (SAS)"=col[14],"mouse (Mmc)"=col[13],"mouse (Mmd)"=col[13],"mouse (Mmm)"=col[13], "mouse (Ms)"=col[3],"bonobo"=col[4],"chimpanzee"=col[5],"orangutan (Sumatran)"=col[8],"orangutan (Bornean)"=col[6],"vaquita"=col[9],"wolf (Broad)"="grey","wolf (UCLA)"="black")

ggplot(niceLabels,aes(y=niceLabel,x="color",color=niceLabel))+
  geom_point(size=5)+
  scale_fill_manual(values=colorList)+
  scale_color_manual(values=colorList)+
  theme_classic()




