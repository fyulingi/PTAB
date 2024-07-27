select ?book ?pages where{
?book <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://yago-knowledge.org/resource/wordnet_book_106410904>.
?book <http://yago-knowledge.org/resource/hasPages> ?pages.
} order by desc(?pages) limit 5