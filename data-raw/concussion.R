## code to prepare `concussion` dataset goes here
library(dplyr)
library(stringr)

ConcData <- read.table(list.files('../data/', 'Symptom', full.names=TRUE),
                       sep='\t', header=TRUE, row.names=1)
colnames(ConcData) <- sub('pcsi_', '', colnames(ConcData))
colnames(ConcData) <- paste0('C', str_extract(colnames(ConcData), '[0-9]+$'),
                             '_', sub('_.*', '', colnames(ConcData)))

head(ConcData)

## Count NAs
apply(ConcData, 1, function(x) sum(is.na(x))) %>% table

## Remove observations with NAs
ConcDataCl <- ConcData[apply(ConcData, 1, function(x) sum(is.na(x))) == 0, ]
class(ConcDataCl)

ConcDataCl <- ConcDataCl + 1
ConcData <- ConcData + 1

allLevs <- paste0('c', 1:7)
concussion <- ConcDataCl %>%
   apply(2, function(x) paste0('c', x)) %>%
   ## apply(2, as.character) %>%
   data.frame(check.names=FALSE) %>%
   mutate_if(is.character, function(x) factor(x, levels = allLevs))
rownames(concussion) <- rownames(ConcDataCl)
head(concussion)
concussion[, 1]

theQuests <- paste0('Q', 1:21, ': ',
                    c('Headache', 'Nausea', 'Balance problems', 'Dizziness',
                      'Fatigue', 'Sleep more', 'Drowsiness',
                      'Sensibility to light', 'Sensibility to noice',
                      'Irritability', 'Sadness', 'Nervousness', 'More emotional',
                      'Feeling slowed down', 'Feeling mentally foggy',
                      'Difficulty concentrating', 'Difficulty remembering',
                      'Visual problem', 'Confusion', 'Feeling clumsy',
                      'Answer slowlier'))
colnames(concussion) <-  theQuests
head(concussion)

usethis::use_data(concussion, overwrite = TRUE)
