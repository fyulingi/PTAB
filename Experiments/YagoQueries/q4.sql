select distinct ?district ?country ?nPeople1 ?nPeople2 ?e ?nPeople3
where{
?district <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://yago-knowledge.org/resource/wordnet_administrative_district_108491826>.
?district <http://yago-knowledge.org/resource/hasNumberOfPeople> ?nPeople1.
?district  <http://yago-knowledge.org/resource/isLocatedIn> ?country.
?country <http://yago-knowledge.org/resource/hasNumberOfPeople> ?nPeople2.
?country <http://yago-knowledge.org/resource/isLocatedIn> ?e.
?e <http://yago-knowledge.org/resource/hasNumberOfPeople> ?nPeople3. 
} order by desc(?nPeople1+?nPeople2+?nPeople3) limit 5