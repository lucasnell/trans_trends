#===load packages=====
library(tidyverse)
library(lubridate)



###################PITFALLS###################################
#====read data====
pit_old <- read_csv("data/myvatn_pitfalls.csv", na = "NA")
pit_new <- read_csv("data/pit_13_19.csv", na = "NA")


all(colnames(pit_old)==colnames(pit_new))

#====find columns that don't match=====
#find columns in old data that don't match new data
colnames(pit_old)[which(!colnames(pit_old) %in% colnames(pit_new))]

colnames(pit_new)[which(!colnames(pit_new) %in% colnames(pit_old))]


#====remove date and daysout====
pit_new <- pit_new %>%
    select(-andate, -daysout)

#check
colnames(pit_new)[which(!colnames(pit_new) %in% colnames(pit_old))]



#====lyco_total====

#check whether they used Lyco_NEW_total or lyco_OLD_total
left_join(pit_old %>% mutate(egsa = as.numeric(egsa)), pit_new) %>%
    filter(Lyco_NEW_total !=lyco_total)
#0 rows

left_join(pit_old %>% mutate(egsa = as.numeric(egsa)), pit_new) %>%
    filter(lyco_OLD_total !=lyco_total)
#66 rows

#thus Lyco_NEW_total = lyco_total

pit_new <- pit_new %>%
    rename(lyco_total = Lyco_NEW_total) %>%
    select(-lyco_OLD_total, -lyco_male, -lyco_female)

#check
colnames(pit_old)[which(!colnames(pit_old) %in% colnames(pit_new))]

colnames(pit_new)[which(!colnames(pit_new) %in% colnames(pit_old))]

#====aran_other====
#confirm aran = aran_other
left_join(pit_old %>% mutate(egsa = as.numeric(egsa)), pit_new) %>%
    filter(aran !=aran_other)

left_join(pit_old %>% mutate(egsa = as.numeric(egsa)), pit_new) %>%
    filter(aran ==aran_other)


pit_new <- pit_new %>%
    rename(aran_other = aran)

#check
colnames(pit_old)[which(!colnames(pit_old) %in% colnames(pit_new))]

colnames(pit_new)[which(!colnames(pit_new) %in% colnames(pit_old))]


#====trom====
left_join(pit_old %>% mutate(egsa = as.numeric(egsa)), pit_new %>% rename(trom2 = trom)) %>%
    filter(trom2 !=trom)

pit_new %>%
    filter(year(coldate)>2017, is.na(trom), !is.na(trom_adult))

pit_new <- pit_new %>%
    mutate(trom = ifelse(trans == "btl" & dist == "5" & coldate == "2018-06-14", 0, trom)) %>%
    select(-trom_adult, -trom_juv)


#check
colnames(pit_old)[which(!colnames(pit_old) %in% colnames(pit_new))]

colnames(pit_new)[which(!colnames(pit_new) %in% colnames(pit_old))]


#====acar_other and acar_total=====
left_join(pit_old %>% mutate(egsa = as.numeric(egsa)), pit_new) %>%
    filter(acar !=acar_other)

pit_old %>%
    filter(acar_total == acar_other + trom)

pit_new <- pit_new %>%
    rename(acar_other = acar) %>%
    mutate(acar_total = acar_other + trom)


#check
colnames(pit_old)[which(!colnames(pit_old) %in% colnames(pit_new))]

colnames(pit_new)[which(!colnames(pit_new) %in% colnames(pit_old))]


#====cole_other====

left_join(pit_old %>% mutate(egsa = as.numeric(egsa)), pit_new) %>%
    filter(cole !=cole_other)


pit_new <- pit_new %>%
    rename(cole_other = cole)

#check
colnames(pit_old)[which(!colnames(pit_old) %in% colnames(pit_new))]

colnames(pit_new)[which(!colnames(pit_new) %in% colnames(pit_old))]


#====final column check====
colnames(pit_old) == colnames(pit_new)

#order to match
pit_new <- pit_new %>%
    select(colnames(pit_old))
#columns now match



#====Find dates that don't match over the overlapping intervals====
pit_old %>% count(year(coldate))
pit_new %>% count(year(coldate))

pit_old %>%
    arrange(coldate, lakeid, trans, dist) %>%
    mutate(year= year(coldate),
           sampleid = paste(trans, dist, coldate)) %>%
    select(year, sampleid)

pit_new %>%
    arrange(coldate, lakeid, trans, dist) %>%
    mutate(year= year(coldate),
           sampleid = paste(trans, dist, coldate)) %>%
    select(year, sampleid)

#remove 2011 from pit_new, because it is not as complete as in pit_old.
#this allows us to focus on the period between 2013 and 2017
pit_new <- pit_new %>%
    filter(year(coldate)!=2011)



pit_old %>%
    filter(year(coldate)>2012) %>%
    count(year(coldate)) %>%
    rename(year = 1, old_n = 2) %>%
    left_join(pit_new %>%
                  filter(year(coldate)<2018) %>%
                  count(year(coldate)) %>%
                  rename(year = 1, new_n = 2))

#2017 shows a different number of rows with more in the new set.


pit_old %>%
    filter(year(coldate)==2017) %>%
    select(coldate) %>%
    unique() %>%
    pull(coldate)

pit_new %>%
    filter(year(coldate)==2017) %>%
    select(coldate) %>%
    unique() %>%
    pull(coldate)
#it looks like the last sampling date is missing in pit_old

olddates <- pit_old %>%
    filter(year(coldate) %in% 2013:2017) %>%
    mutate(sampleid = paste(trans, dist, setdate)) %>%
    pull(sampleid)

newdates <- pit_new %>%
    filter(year(coldate) %in% 2013:2017) %>%
    mutate(sampleid = paste(trans, dist, setdate)) %>%
    pull(sampleid)

which(!olddates %in% newdates)
newdates[which(!newdates %in% olddates)]
#again all match but the one day.

olddates <- pit_old %>%
    filter(year(coldate) %in% 2013:2017) %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    pull(sampleid)

newdates <- pit_new %>%
    filter(year(coldate) %in% 2013:2017) %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    pull(sampleid)

which(!olddates %in% newdates)
newdates[which(!newdates %in% olddates)]

#just the one coldate



# Now that we have confirmed the sampleids match, check that the data match

#==== check arthropod data ====

#in pit old, there were some days where egsa was listed as "x" I'm going to make these NA
pit_old %>%
    filter(!is.na(egsa) & is.na(as.numeric(egsa))) %>%
    select(trans, dist, coldate, egsa)


pit_old <- pit_old %>%
    mutate(egsa = as.numeric(egsa))

pitcomp1 <- pit_new %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    select(sampleid)

pitcomp2 <- pit_old %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    select(sampleid)

pitcomp1$new_total <- rowSums(pit_new[,6:33], na.rm = TRUE)

pitcomp2$old_total <- rowSums(pit_old[,6:33], na.rm = TRUE)


pitcomp <- pitcomp1 %>% left_join(pitcomp2)

pitcomp %>% filter(!str_detect(sampleid, "2018"), !str_detect(sampleid, "2019")) %>% print(n=1000)

#totals all match so things look good.


#the comments seem to not match
pit_new %>%
    filter(year(coldate) %in% 2013:2017) %>%
    rename(new_comments = comments,
           new_comments2 = `comments 2`) %>%
    full_join(pit_old %>%
                  filter(year(coldate) %in% 2013:2017)) %>%
    filter(comments!=new_comments) %>%
    select(trans, dist, coldate, comments, new_comments) %>%
    print(n = 1000)


pit_new %>%
    filter(year(coldate) %in% 2013:2017) %>%
    rename(new_comments = comments,
           new_comments2 = `comments 2`) %>%
    full_join(pit_old %>%
                  filter(year(coldate) %in% 2013:2017)) %>%
    filter(`comments 2`!=new_comments2) %>%
    select(trans, dist, coldate, `comments 2`, new_comments2) %>%
    print(n = 1000)

#check ones that says sample lost
pit_old %>% filter(coldate == "2017-06-14", trans == "btl", dist == 150) %>% View
pit_old %>% filter(coldate == "2016-06-17", trans == "hag", dist == 150) %>% View

#there appears to be some sort of frame shift error in the old pitfall data comments
#I will join using the new 2013 to 2019 data.


#====Join Pitfall data====
pit_full <- pit_old %>%
    filter(year(coldate)<2013) %>%
    full_join(pit_new) %>%
    arrange(coldate, lakeid, trans, dist)



pit_full %>%
    gather(var, val, aran_other:clmb) %>%
    ggplot(aes(y = val, x = yday(coldate), col = year(coldate)))+
    geom_line(aes(group = interaction(trans, dist)))+
    facet_wrap(~var)

###################PITFALLS###################################
#====read data====
inf_old <- read_csv("data/myvatn_infalls.csv", na = c("na", "NA"), col_types = "ccdDDddddddddddcc")
inf_new <- read_csv("data/inf_13_19.csv",  na = c("na", "NA"), col_types = "ccdDDDddddddddddcc")


#====check match in columns====
unique(inf_new$X18)
inf_new <- inf_new %>%
    select(-X18)
all(colnames(inf_old)==colnames(inf_new))

colnames(inf_old)[which(!colnames(inf_old) %in% colnames(inf_new))]
colnames(inf_new)[which(!colnames(inf_new) %in% colnames(inf_old))]

#===remove andate and add comments.2====
inf_new <- inf_new %>%
    select(-andate) %>%
    mutate(comments.2 = "")

#check
colnames(inf_old)[which(!colnames(inf_old) %in% colnames(inf_new))]
colnames(inf_new)[which(!colnames(inf_new) %in% colnames(inf_old))]


#====check dates====
inf_old %>%
    count(year(coldate)) %>%
    rename(year = 1, old_n =2) %>%
    full_join(inf_new %>%
                  count(year(coldate)) %>%
                  rename(year = 1, new_n = 2))

olddates <- inf_old %>%
    filter(year(coldate) %in% 2013:2017) %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    select(sampleid) %>%
    unique() %>%
    pull(sampleid)

newdates <- inf_new %>%
    filter(year(coldate) %in% 2013:2017) %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    select(sampleid) %>%
    unique() %>%
    pull(sampleid)




olddates[which(!olddates %in% newdates)]
newdates[which(!newdates %in% olddates)]


inf_old %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    filter(sampleid %in% olddates[which(!olddates %in% newdates)]) %>%
    as.data.frame()


inf_new %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    filter(sampleid %in% newdates[which(!newdates %in% olddates)])


inf_new %>% filter(trans == "hel2", dist == "500", year(coldate)==2017)



#fix days from 2017 with reversed set and collect dates
inf_new <- inf_new %>%
    mutate(setdate = if_else(trans == "hel2" & setdate == "2017-06-29" & coldate == "2017-06-14", as.Date("2017-06-14"), setdate),
           coldate = if_else(trans == "hel2" & setdate == "2017-06-14" & coldate == "2017-06-14", as.Date("2017-06-29"), coldate))



#check
olddates <- inf_old %>%
    filter(year(coldate) %in% 2013:2017) %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    select(sampleid) %>%
    unique() %>%
    pull(sampleid)

newdates <- inf_new %>%
    filter(year(coldate) %in% 2013:2017) %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    select(sampleid) %>%
    unique() %>%
    pull(sampleid)

olddates[which(!olddates %in% newdates)]
newdates[which(!newdates %in% olddates)]


inf_old %>%
    filter(year(coldate) %in% 2013:2017) %>%
    count(coldate) %>%
    rename(old_n = 2) %>%
    full_join(inf_new %>%
                  filter(year(coldate)%in% 2013:2017) %>%
                  count(coldate) %>%
                  rename(new_n = 2))

full_join(inf_old, inf_new)



#remove duplicated and incorrect day
inf_new <- inf_new %>%
    filter(!(trans == "btl" & dist == "5" & coldate == "2013-08-14" & setdate == "2013-06-05"))



#===check row data===


infcomp1 <- inf_new %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    select(sampleid)

infcomp2 <- inf_old %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    select(sampleid)

infcomp1$new_total <- inf_new$bgch + inf_new$smch + inf_new$simu

infcomp2$old_total <- inf_old$bgch + inf_old$smch + inf_old$simu


infcomp <- infcomp1 %>% left_join(infcomp2)

infcomp %>% filter(new_total!=old_total)

infcomp %>%
    filter(!str_detect(sampleid, "2018"), !str_detect(sampleid, "2019")) %>%
    print(n=1000)




infcomp1$fract_new <- inf_new$fract_bgch + inf_new$fract_smch + inf_new$fract_simu

infcomp2$fract_old <- inf_old$fract_bgch + inf_old$fract_smch + inf_old$fract_simu


infcomp <- infcomp1 %>% left_join(infcomp2)

infcomp %>% filter(new_total!=old_total)

infcomp %>%
    filter(!str_detect(sampleid, "2018"), !str_detect(sampleid, "2019")) %>%
    print(n=1000)
#totals all match so things look good.

#check comments
inf_new %>%
    filter(year(coldate)%in% 2013:2017,
           coldate!="2017-08-18") %>%
    full_join(inf_old %>% filter(year(coldate) %in% 2013:2017) %>%
                  rename(commentsold = comments, comments2old = comments.2)) %>%
    arrange(coldate, lakeid, trans, dist) %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    count(sampleid) %>%
    arrange(desc(n))

#there are some duplicated rows

inf_new %>% filter(trans == "hel1", dist == "500", coldate == "2017-06-14") %>%
    as.data.frame()

inf_old %>% filter(trans == "hel1", dist == "500", coldate == "2017-06-14") %>%
    as.data.frame()

inf_new <- unique(inf_new)
inf_old <- unique(inf_old)


inf_new %>%
    filter(year(coldate)%in% 2013:2017,
           coldate!="2017-08-18") %>%
    full_join(inf_old %>% filter(year(coldate) %in% 2013:2017) %>%
                  rename(commentsold = comments, comments2old = comments.2)) %>%
    arrange(coldate, lakeid, trans, dist) %>%
    mutate(sampleid = paste(trans, dist, coldate)) %>%
    count(sampleid) %>%
    unique() %>%
    arrange(desc(n))

inf_new %>% filter(trans == "non", dist == "50", coldate == "2015-07-30") %>%
    as.data.frame()

inf_old %>% filter(trans == "non", dist == "50", coldate == "2015-07-30") %>%
    as.data.frame()



inf_full <- inf_old %>%
    filter(year(coldate)<2013) %>%
    full_join(inf_new) %>%
    arrange(coldate, lakeid, trans, dist)





######Final checks####################
pit <- pit_full
inf <- inf_full


pit %>%
    group_by(trans) %>%
    ggplot(aes(coldate, fill = trans)) +
    facet_wrap(~year(coldate), scales = "free_x")+
    geom_histogram()

pit %>%
    group_by(trans) %>%
    ggplot(aes(coldate, fill = trans)) +
    facet_wrap(~year(coldate), scales = "free_x")+
    geom_histogram()

pit %>%
    filter(trans == "btl", dist == 200)

pit %>% filter(coldate == "2012-07-03")


d <- crossing(trans = unique(pit$trans), dist = unique(pit$dist)) %>%
    mutate(levels = paste(trans, dist))

d <- pit %>%
    filter(!is.na(lakeid), trans !="kal2") %>%
    select(lakeid, trans, dist) %>%
    unique()%>%
    mutate(levels = paste(trans, dist)) %>%
    filter(levels!= "btl 200") %>%
    arrange(lakeid, trans, dist)

pitplot <- pit %>%
    filter(trans!="kal2") %>%
    mutate(trap = paste(trans, dist),
           trap = factor(trap, levels = d$levels)) %>%
    filter(trap!="btl 200")

pitplot$year = year(pitplot$coldate)
year(pitplot$coldate) <- 3000
year(pitplot$setdate) <- 3000

#pit lollipop plot
pitplot%>%
    ggplot(aes(x = trap)) +
    facet_wrap(~year)+
    geom_segment(aes(y = setdate, yend = coldate, xend =trap, color = trans))+
    geom_point(aes(y = coldate, color = trans))+
    scale_x_discrete(labels = d$dist)+
    theme_bw()+
    scale_y_date(date_breaks = c("1 month"), date_labels = "%b")+
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_blank())

#looks good

inf %>%
    group_by(trans) %>%
    ggplot(aes(coldate, fill = trans)) +
    facet_wrap(~year(coldate), scales = "free_x")+
    geom_histogram()

inf %>%
    group_by(trans) %>%
    ggplot(aes(coldate, fill = trans)) +
    facet_wrap(~year(coldate), scales = "free_x")+
    geom_histogram()

e <- crossing(trans = unique(inf$trans), dist = unique(inf$dist)) %>%
    mutate(levels = paste(trans, dist))

e <- inf %>%
    filter(!is.na(lakeid), trans !="kal2") %>%
    select(lakeid, trans, dist) %>%
    unique()%>%
    mutate(levels = paste(trans, dist)) %>%
    filter(levels!= "btl 200") %>%
    arrange(lakeid, trans, dist)

infplot <- inf %>%
    filter(trans!="kal2") %>%
    mutate(trap = paste(trans, dist),
           trap = factor(trap, levels = e$levels)) %>%
    filter(trap!="btl 200")

infplot$year = year(infplot$coldate)
year(infplot$coldate) <- 3000
year(infplot$setdate) <- 3000

infplot%>%
    ggplot(aes(x = trap)) +
    facet_wrap(~year)+
    geom_segment(aes(y = setdate, yend = coldate, xend =trap, color = trans))+
    geom_point(aes(y = coldate, color = trans))+
    scale_x_discrete(labels = e$dist)+
    theme_bw()+
    scale_y_date(date_breaks = c("1 month"), date_labels = "%b")+
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_blank())



inf %>%
    filter(trans == "skf", year(coldate) == 2012, month(setdate)==6) %>%
    as.data.frame()
#I suspect that the distances should be 20, 50, 200, and 350 rather than 5, 50, 150, 500
#even though previous years the distances are 5, 50, 150, 500

#I will change those and leave comments

inf <- inf %>%
    mutate(comments.2 =ifelse(trans == "skf" & coldate == "2012-06-19", "JB: These were originally labeled as 5, 50, 150, 500 changed to 20, 50, 200, and 350 to match other sampling events in 2012. 8 Oct 2021", comments.2),
           dist = ifelse(trans == "skf" & coldate == "2012-06-19" & dist == 5, 20, dist),
           dist = ifelse(trans == "skf" & coldate == "2012-06-19" & dist == 150, 200, dist),
           dist = ifelse(trans == "skf" & coldate == "2012-06-19" & dist == 500, 350, dist))


#====There are some lingering issues with duplicates on dates====


inf %>%
    unique() %>%
    group_by(trans) %>%
    count(dist, setdate, coldate) %>%
    filter(n>1) %>%
    select(-n)

inf %>%
    filter(trans == "hag", dist == 5, year(coldate) == 2012) %>%
    as.data.frame()

inf %>%
    filter(trans == "hel2", dist == 50, year(coldate) == 2018) %>%
    as.data.frame()
#these differ only in comments

#this is
inf %>%
    filter(trans == "non",  year(coldate) == 2012) %>%
    as.data.frame()


inf %>%
    filter(trans == "kal", year(coldate) == 2012) %>%
    as.data.frame()

inf %>%
    filter(trans == "non", dist == 150, year(coldate) == 2008) %>%
    as.data.frame()


###################WRITE CSVS####################
# write_csv(pit, "data/myvatn_pitfalls_08-19.csv")
# write_csv(inf, "data/myvatn_infalls_08-19.csv")
