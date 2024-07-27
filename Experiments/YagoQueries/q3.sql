select ?politician ?place ?area ?gdp
 where{
?politician <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://yago-knowledge.org/resource/wordnet_politician_110451263>.
?politician <http://yago-knowledge.org/resource/hasWikipediaArticleLength> ?len1.
?politician <http://yago-knowledge.org/resource/wasBornIn> ?place.
?place <http://yago-knowledge.org/resource/hasWikipediaArticleLength> ?len2.
?place <http://yago-knowledge.org/resource/isLocatedIn> ?area.
?area <http://yago-knowledge.org/resource/hasGDP> ?gdp.
}order by ?gdp limit 5