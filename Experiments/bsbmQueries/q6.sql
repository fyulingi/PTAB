select ?v0 ?v1 ?v3 ?v4 ?v5 ?v6 ?N1 ?N2 ?N3 ?N4 where
{
	?v0 <http://www.w3.org/2000/01/rdf-schema#subClassOf> ?v1 .
	?v3 <http://purl.org/dc/elements/1.1/publisher> ?v4 .
	?v3 <http://www.w3.org/2000/01/rdf-schema#subClassOf> ?v1 .
	?v5 <http://purl.org/dc/elements/1.1/publisher> ?v6 .
	?pro <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?v3.
	?pro <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric1> ?N1.
	?pro <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric2> ?N2.
	?pro <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric3> ?N3.
	?pro <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric4> ?N4.
}order by (?N1+?N2+?N3+?N4) limit 5

