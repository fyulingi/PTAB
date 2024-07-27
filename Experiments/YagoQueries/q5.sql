select * where{
?war <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://yago-knowledge.org/resource/wikicategory_1870s_conflicts>.
?country <http://yago-knowledge.org/resource/participatedIn> ?war.
?area <http://yago-knowledge.org/resource/isLocatedIn> ?country.
?district <http://yago-knowledge.org/resource/isLocatedIn> ?area.
?district <http://yago-knowledge.org/resource/hasGeonamesEntityId> ?id.
} order by ?id limit 5