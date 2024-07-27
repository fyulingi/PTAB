select distinct ?district ?country ?nPeople1 ?nPeople2 ?e1 ?nPeople3 ?e2 ?nPeople4
where{
?district <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://yago-knowledge.org/resource/wordnet_administrative_district_108491826>.
?district <http://yago-knowledge.org/resource/hasNumberOfPeople> ?nPeople1.
?district  <http://yago-knowledge.org/resource/isLocatedIn> ?country.
?country <http://yago-knowledge.org/resource/hasNumberOfPeople> ?nPeople2.
?country <http://yago-knowledge.org/resource/isLocatedIn> ?e1.
?country <http://yago-knowledge.org/resource/isLocatedIn> ?e2.
?e1 <http://yago-knowledge.org/resource/hasNumberOfPeople> ?nPeople3. 
?e2 <http://yago-knowledge.org/resource/hasNumberOfPeople> ?nPeople4.
?e1 <http://yago-knowledge.org/resource/hasWebsite> ?website.
?e2 <http://yago-knowledge.org/resource/exports>  <http://yago-knowledge.org/resource/wordnet_chemical_114806838>.
} order by desc(?nPeople1+?nPeople2+?nPeople3+?nPeople4) 
limit 5