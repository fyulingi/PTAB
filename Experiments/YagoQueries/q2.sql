select ?politician ?palce ?len1 ?len2  where{
?politician <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://yago-knowledge.org/resource/wordnet_politician_110451263>.
?politician <http://yago-knowledge.org/resource/hasWikipediaArticleLength> ?len1.
?politician <http://yago-knowledge.org/resource/wasBornIn> ?place.
?place <http://yago-knowledge.org/resource/hasWikipediaArticleLength> ?len2.
} order by desc(?len1 + ?len2)
limit 5