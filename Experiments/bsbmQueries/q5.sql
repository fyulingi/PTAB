select ?v0 ?v1 ?v2 ?v3 ?v4 ?v5 ?v6 ?v7 ?v8 
?instances0 ?instances5 
?N0 ?N5 
where
{
	?v0 <http://purl.org/dc/elements/1.1/date> "2000-07-17"^^<http://www.w3.org/2001/XMLSchema#date> .
	?v0 <http://purl.org/dc/elements/1.1/publisher> ?v1 .
	?v0 <http://www.w3.org/2000/01/rdf-schema#label> ?v2 .
	?v0 <http://www.w3.org/2000/01/rdf-schema#comment> ?v8 .
	?v3 <http://www.w3.org/2000/01/rdf-schema#subClassOf> <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/instances/ProductType2> .
	?v3 <http://purl.org/dc/elements/1.1/publisher> ?v1 .
	?v3 <http://www.w3.org/2000/01/rdf-schema#label> ?v4 .
	?v3 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?v7 .
	?v5 <http://purl.org/dc/elements/1.1/publisher> ?v6 .
	?v5 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?v7 .
	?instances5 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?v5.		
	?instances5 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric1> ?N5.
}
order by desc(?N5)
limit 5