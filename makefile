#help for make
#http://www.cnblogs.com/wang_yb/p/3990952.html
#https://segmentfault.com/a/1190000000349917
#http://blog.csdn.net/cuiyifang/article/details/7910268

#to use gprof to analyse efficience of the program:
#http://blog.chinaunix.net/uid-25194149-id-3215487.html

#to use gcov and lcov
#Notice that optimization should not be used here
#http://blog.163.com/bobile45@126/blog/static/96061992201382025729313/
#gcov -a main.cpp
#lcov --directory .   --capture --output-file dig.info 
#genhtml --output-directory . --frames --show-details dig.info 

#to use doxygen+graphviz+htmlhelp to generate document from source code:
#http://www.doxygen.nl/
#(also include good comments norm)
#http://blog.csdn.net/u010740725/article/details/51387810

#CC=$(shell which clang 2>/dev/null || which gcc)
#ccache, readline, gcov lcov
#http://blog.csdn.net/u012421852/article/details/52138960
#
# How to speed up the compilation
# https://blog.csdn.net/a_little_a_day/article/details/78251928
# use make -j4, if error then use make utilizing only one thread
#use -j8 or higher may cause error
#http://blog.csdn.net/cscrazybing/article/details/50789482
#http://blog.163.com/liuhonggaono1@126/blog/static/10497901201210254622141/

#compile parameters

# WARN: maybe difficult to install ccache in some systems
#CC = ccache g++
CXX = g++

#the optimazition level of gcc/g++
#http://blog.csdn.net/hit_090420216/article/details/44900215
#NOTICE: -O2 is recommended, while -O3(add loop-unroll and inline-function) is dangerous
#when developing, not use -O because it will disturb the normal 
#routine. use it for test and release.
CFLAGS = -c -Wall -O2 -pthread -std=c++11 -Werror=return-type
EXEFLAG = -O2 -pthread -std=c++11 -Werror=return-type
#-coverage for debugging
#CFLAGS = -c -Wall -pthread -g3 -std=c++11  -gdwarf-2 -Werror=return-type
#EXEFLAG = -pthread -g3 -std=c++11 -gdwarf-2 -Werror=return-type
#-coverage for debugging and with performance
# CFLAGS = -c -Wall -pthread -g3 -std=c++11  -gdwarf-2 -pg
# EXEFLAG = -pthread -g3 -std=c++11 -gdwarf-2 -pg

#add -lreadline [-ltermcap] if using readline or objs contain readline
# library = -lreadline -L./lib -L/usr/local/lib -lantlr -lgcov -lboost_thread -lboost_filesystem -lboost_system -lboost_regex -lpthread -I/usr/local/include/boost -lcurl
library = -lreadline -L./lib -L/usr/local/lib -L/usr/lib/ -lantlr4-runtime -lgcov -lboost_thread -lboost_filesystem -lboost_system -lboost_regex -lpthread -I/usr/local/include/boost -lcurl
#used for parallelsort
openmp = -fopenmp -march=native
# library = -ltermcap -lreadline -L./lib -lantlr -lgcov
def64IO = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE

# paths

objdir = .objs/

exedir = bin/

testdir = scripts/

lib_antlr = lib/libantlr4-runtime.a

api_cpp = api/socket/cpp/lib/libgstoreconnector.a

api_java = api/socket/java/lib/GstoreJavaAPI.jar

# objects

#sstreeobj = $(objdir)Tree.o $(objdir)Storage.o $(objdir)Node.o $(objdir)IntlNode.o $(objdir)LeafNode.o $(objdir)Heap.o 
sitreeobj = $(objdir)SITree.o $(objdir)SIStorage.o $(objdir)SINode.o $(objdir)SIIntlNode.o $(objdir)SILeafNode.o $(objdir)SIHeap.o 
istreeobj = $(objdir)ISTree.o $(objdir)ISStorage.o $(objdir)ISNode.o $(objdir)ISIntlNode.o $(objdir)ISLeafNode.o $(objdir)ISHeap.o 
ivtreeobj = $(objdir)IVTree.o $(objdir)IVStorage.o $(objdir)IVNode.o $(objdir)IVIntlNode.o $(objdir)IVLeafNode.o $(objdir)IVHeap.o 
ivarrayobj = $(objdir)IVArray.o $(objdir)IVEntry.o $(objdir)IVBlockManager.o
isarrayobj = $(objdir)ISArray.o $(objdir)ISEntry.o $(objdir)ISBlockManager.o

kvstoreobj = $(objdir)KVstore.o $(sitreeobj) $(istreeobj) $(ivtreeobj) $(ivarrayobj) $(isarrayobj) #$(sstreeobj)

utilobj = $(objdir)Util.o $(objdir)Bstr.o $(objdir)Stream.o $(objdir)Triple.o $(objdir)BloomFilter.o $(objdir)VList.o \
			$(objdir)EvalMultitypeValue.o $(objdir)IDTriple.o $(objdir)Version.o $(objdir)Transaction.o $(objdir)Latch.o $(objdir)IPWhiteList.o \
			$(objdir)IPBlackList.o $(objdir)OrderedVector.o

topkobj = $(objdir)DynamicTrie.o $(objdir)CompressedVector.o $(objdir)DPBList.o $(objdir)DPPList.o $(objdir)HopIndex.o $(objdir)Pool.o $(objdir)TopKUtil.o $(objdir)DPBTopKUtil.o \
			$(objdir)DPPTopKUtil.o $(objdir)TopKSearchPlan.o  $(objdir)GlobalQueue.o $(objdir)DPPFQIterator.o \
			$(objdir)DPPFRIterator.o $(objdir)DPPOWIterator.o $(objdir)Subspace.o $(objdir)kTPMUtil.o \
			$(objdir)SAEUtil.o $(objdir)SpaceSAE.o $(objdir)DfsFQIterator.o\
			$(objdir)DfsFQIteratorOld.o $(objdir)DfsFRIterator.o $(objdir)DfsOWIterator.o $(objdir)DfsList.o $(objdir)DfsUtil.o \
			$(objdir)DfsUtilOld.o $(objdir)EagerUtil.o $(objdir)Take2Util.o $(objdir)TakeAllUtil.o \
			$(objdir)EagerSubspace.o $(objdir)Take2Subspace.o $(objdir)TakeAllSubspace.o

queryobj = $(objdir)SPARQLquery.o $(objdir)BasicQuery.o $(objdir)ResultSet.o  $(objdir)IDList.o $(objdir)QueryPlan.o\
		   $(objdir)Varset.o $(objdir)QueryTree.o $(objdir)TempResult.o $(objdir)QueryCache.o $(objdir)GeneralEvaluation.o \
		   $(objdir)PathQueryHandler.o $(objdir)BGPQuery.o $(topkobj)

signatureobj = $(objdir)SigEntry.o $(objdir)Signature.o

vstreeobj = $(objdir)VSTree.o $(objdir)EntryBuffer.o $(objdir)LRUCache.o $(objdir)VNode.o

stringindexobj = $(objdir)StringIndex.o

parserobj = $(objdir)RDFParser.o $(objdir)SPARQLParser.o $(objdir)DBparser.o \
			$(objdir)SPARQLLexer.o $(objdir)TurtleParser.o $(objdir)QueryParser.o

serverobj = $(objdir)Operation.o $(objdir)Server.o $(objdir)Client.o $(objdir)Socket.o

# httpobj = $(objdir)client_http.hpp.gch $(objdir)server_http.hpp.gch

databaseobj = $(objdir)Statistics.o $(objdir)Database.o $(objdir)Join.o $(objdir)Strategy.o \
 $(objdir)CSR.o $(objdir)Txn_manager.o $(objdir)TableOperator.o $(objdir)PlanTree.o  \
 $(objdir)PlanGenerator.o $(objdir)Executor.o $(objdir)Optimizer.o $(objdir)ResultTrigger.o

trieobj = $(objdir)Trie.o $(objdir)TrieNode.o

objfile = $(kvstoreobj) $(vstreeobj) $(stringindexobj) $(parserobj) $(serverobj) $(httpobj) $(databaseobj) \
		  $(utilobj) $(signatureobj) $(queryobj) $(trieobj)

	 
inc = -I./tools/antlr4-cpp-runtime-4/runtime/src
#-I./usr/local/include/boost/


#auto generate dependencies
# http://blog.csdn.net/gmpy_tiger/article/details/51849474
# http://blog.csdn.net/jeffrey0000/article/details/12421317

#gtest

TARGET = $(exedir)gexport $(exedir)gbuild $(exedir)gserver $(exedir)gserver_backup_scheduler $(exedir)gclient $(exedir)gquery \
 		 $(exedir)gconsole $(api_java) $(exedir)gadd $(exedir)gsub $(exedir)ghttp $(exedir)gmonitor $(exedir)gshow $(exedir)shutdown \
 		 $(exedir)indexBuilder $(exedir)ginit $(exedir)gdrop $(testdir)update_test $(testdir)dataset_test $(testdir)transaction_test \
 		 $(testdir)run_transaction $(exedir)gbackup $(exedir)grestore $(exedir)gpara $(exedir)rollback \
 		 $(exedir)TestTotalTime $(exedir)TestKValueOneQuery $(exedir)TestRootNode $(exedir)TestSearchSpace $(exedir)TestCycle

all: $(TARGET)
	@echo "Compilation ends successfully!"
	@bash scripts/init.sh

#BETTER: use for loop to reduce the lines
#NOTICE: g++ -MM will run error if linking failed, like Database.h/../SparlParser.h/../antlr3.h

#executables begin

#NOTICE:not include g*.o in objfile due to multiple definitions of main()

$(exedir)gexport: $(lib_antlr) $(objdir)gexport.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)gexport $(objdir)gexport.o $(objfile) $(library) $(openmp)

$(exedir)gdrop: $(lib_antlr) $(objdir)gdrop.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)gdrop $(objdir)gdrop.o $(objfile) $(library) $(openmp)

$(exedir)ginit: $(lib_antlr) $(objdir)ginit.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)ginit $(objdir)ginit.o $(objfile) $(library) $(openmp)

$(exedir)shutdown: $(lib_antlr) $(objdir)shutdown.o $(objfile) $(api_cpp)
	$(CXX) $(EXEFLAG) -o $(exedir)shutdown $(objdir)shutdown.o $(objfile) $(openmp) -L./api/http/cpp/lib -lclient $(library)

$(exedir)gmonitor: $(lib_antlr) $(objdir)gmonitor.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)gmonitor $(objdir)gmonitor.o $(objfile) $(library) $(openmp)

$(exedir)gshow: $(lib_antlr) $(objdir)gshow.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)gshow $(objdir)gshow.o $(objfile) $(library) $(openmp)

$(exedir)gbuild: $(lib_antlr) $(objdir)gbuild.o $(objfile) 
	$(CXX) $(EXEFLAG) -o $(exedir)gbuild $(objdir)gbuild.o $(objfile) $(library) $(openmp)

$(exedir)gquery: $(lib_antlr) $(objdir)gquery.o $(objfile) 
	$(CXX) $(EXEFLAG) -o $(exedir)gquery $(objdir)gquery.o $(objfile) $(library) $(openmp)

$(exedir)gserver: $(lib_antlr) $(objdir)gserver.o $(objfile) 
	$(CXX) $(EXEFLAG) -o $(exedir)gserver $(objdir)gserver.o $(objfile) $(library) $(openmp)

$(exedir)gserver_backup_scheduler: $(lib_antlr) $(objdir)gserver_backup_scheduler.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)gserver_backup_scheduler $(objdir)gserver_backup_scheduler.o $(objfile) $(library) $(openmp)

$(exedir)gclient: $(lib_antlr) $(objdir)gclient.o $(objfile) 
	$(CXX) $(EXEFLAG) -o $(exedir)gclient $(objdir)gclient.o $(objfile) $(library) $(openmp)

$(exedir)gconsole: $(lib_antlr) $(objdir)gconsole.o $(objfile) $(api_cpp)
	$(CXX) $(EXEFLAG) -o $(exedir)gconsole $(objdir)gconsole.o $(objfile) $(library) -L./api/socket/cpp/lib -lgstoreconnector $(openmp)

$(exedir)ghttp: $(lib_antlr) $(objdir)ghttp.o ./Server/server_http.hpp ./Server/client_http.hpp $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)ghttp $(objdir)ghttp.o $(objfile) $(library) $(inc) -DUSE_BOOST_REGEX $(openmp)

$(exedir)gbackup: $(lib_antlr) $(objdir)gbackup.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)gbackup $(objdir)gbackup.o $(objfile) $(library) $(openmp)

$(exedir)grestore: $(lib_antlr) $(objdir)grestore.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)grestore $(objdir)grestore.o $(objfile) $(library) $(openmp)

$(exedir)gpara: $(lib_antlr) $(objdir)gpara.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)gpara $(objdir)gpara.o $(objfile) $(library) $(openmp) -L./api/http/cpp/lib -lclient $(library)

$(exedir)rollback: $(lib_antlr) $(objdir)rollback.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)rollback $(objdir)rollback.o $(objfile) $(library) $(openmp) -L./api/http/cpp/lib -lclient $(library)

$(testdir)update_test: $(lib_antlr) $(objdir)update_test.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(testdir)update_test $(objdir)update_test.o $(objfile) $(library) $(openmp)

$(testdir)dataset_test: $(lib_antlr) $(objdir)dataset_test.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(testdir)dataset_test $(objdir)dataset_test.o $(objfile) $(library) $(openmp)

$(testdir)transaction_test: $(lib_antlr) $(objdir)transaction_test.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(testdir)transaction_test $(objdir)transaction_test.o $(objfile) $(library) $(openmp)

$(testdir)run_transaction: $(lib_antlr) $(objdir)run_transaction.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(testdir)run_transaction $(objdir)run_transaction.o $(objfile) $(library) $(openmp) -L./api/http/cpp/lib -lclient $(library)


$(exedir)TestTotalTime: $(lib_antlr) $(objdir)TestTotalTime.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)TestTotalTime $(objdir)TestTotalTime.o $(objfile) $(library) $(openmp)

$(exedir)TestKValueOneQuery: $(lib_antlr) $(objdir)TestKValueOneQuery.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)TestKValueOneQuery $(objdir)TestKValueOneQuery.o $(objfile) $(library) $(openmp)

$(exedir)TestRootNode: $(lib_antlr) $(objdir)TestRootNode.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)TestRootNode $(objdir)TestRootNode.o $(objfile) $(library) $(openmp)

$(exedir)TestSearchSpace: $(lib_antlr) $(objdir)TestSearchSpace.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)TestSearchSpace $(objdir)TestSearchSpace.o $(objfile) $(library) $(openmp)

$(exedir)TestCycle: $(lib_antlr) $(objdir)TestCycle.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)TestCycle $(objdir)TestCycle.o $(objfile) $(library) $(openmp)
#executables end


#objects in Main/ begin

$(objdir)gexport.o: Main/gexport.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gexport.cpp $(inc) -o $(objdir)gexport.o $(openmp)

$(objdir)gdrop.o: Main/gdrop.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gdrop.cpp $(inc) -o $(objdir)gdrop.o $(openmp)

$(objdir)ginit.o: Main/ginit.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/ginit.cpp $(inc) -o $(objdir)ginit.o $(openmp)

$(objdir)shutdown.o: Main/shutdown.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS)	Main/shutdown.cpp $(inc) -o $(objdir)shutdown.o $(openmp)

$(objdir)gmonitor.o: Main/gmonitor.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gmonitor.cpp $(inc) -o $(objdir)gmonitor.o $(openmp)

$(objdir)gshow.o: Main/gshow.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gshow.cpp $(inc) -o $(objdir)gshow.o $(openmp)

$(objdir)gbuild.o: Main/gbuild.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gbuild.cpp $(inc) -o $(objdir)gbuild.o $(openmp)
	
$(objdir)gquery.o: Main/gquery.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gquery.cpp $(inc) -o $(objdir)gquery.o $(openmp) #-DREADLINE_ON
	#add -DREADLINE_ON if using readline

$(objdir)gserver.o: Main/gserver.cpp Server/Server.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gserver.cpp $(inc) -o $(objdir)gserver.o $(openmp)

$(objdir)gserver_backup_scheduler.o: Main/gserver_backup_scheduler.cpp Server/Server.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gserver_backup_scheduler.cpp $(inc) -o $(objdir)gserver_backup_scheduler.o $(openmp)

$(objdir)gclient.o: Main/gclient.cpp Server/Client.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gclient.cpp $(inc) -o $(objdir)gclient.o $(openmp) #-DREADLINE_ON

$(objdir)gconsole.o: Main/gconsole.cpp Database/Database.h Util/Util.h api/socket/cpp/src/GstoreConnector.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gconsole.cpp $(inc) -o $(objdir)gconsole.o -I./api/socket/cpp/src/ $(openmp) #-DREADLINE_ON

$(objdir)ghttp.o: Main/ghttp.cpp Server/server_http.hpp Server/client_http.hpp Database/Database.h Database/Txn_manager.h Util/Util.h Util/IPWhiteList.h Util/IPBlackList.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/ghttp.cpp $(inc) -o $(objdir)ghttp.o -DUSE_BOOST_REGEX $(def64IO) $(openmp)

$(objdir)gbackup.o: Main/gbackup.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gbackup.cpp $(inc) -o $(objdir)gbackup.o $(openmp)

$(objdir)grestore.o: Main/grestore.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/grestore.cpp $(inc) -o $(objdir)grestore.o $(openmp)

$(objdir)gpara.o: Main/gpara.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/gpara.cpp $(inc) -o $(objdir)gpara.o $(openmp)

$(objdir)rollback.o: Main/rollback.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Main/rollback.cpp $(inc) -o $(objdir)rollback.o $(openmp)

$(objdir)TestTotalTime.o: Experiments/TestTotalTime.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Experiments/TestTotalTime.cpp $(inc) -o $(objdir)TestTotalTime.o $(openmp)

$(objdir)TestKValueOneQuery.o: Experiments/TestKValueOneQuery.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Experiments/TestKValueOneQuery.cpp $(inc) -o $(objdir)TestKValueOneQuery.o $(openmp)

$(objdir)TestRootNode.o: Experiments/TestRootNode.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Experiments/TestRootNode.cpp $(inc) -o $(objdir)TestRootNode.o $(openmp)

$(objdir)TestSearchSpace.o: Experiments/TestSearchSpace.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Experiments/TestSearchSpace.cpp $(inc) -o $(objdir)TestSearchSpace.o $(openmp)

$(objdir)TestCycle.o: Experiments/TestCycle.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) Experiments/TestCycle.cpp $(inc) -o $(objdir)TestCycle.o $(openmp)

#objects in Main/ end

#objects in scripts/ begin

$(objdir)update_test.o: $(testdir)update_test.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) $(testdir)update_test.cpp $(inc) -o $(objdir)update_test.o $(openmp)

$(objdir)dataset_test.o: $(testdir)dataset_test.cpp Database/Database.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) $(testdir)dataset_test.cpp $(inc) -o $(objdir)dataset_test.o $(openmp)

$(objdir)transaction_test.o: $(testdir)transaction_test.cpp Database/Database.h Database/Txn_manager.h Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) $(testdir)transaction_test.cpp $(inc) -o $(objdir)transaction_test.o $(openmp)

$(objdir)run_transaction.o: $(testdir)run_transaction.cpp Util/Util.h $(lib_antlr)
	$(CXX) $(CFLAGS) $(testdir)run_transaction.cpp $(inc) -o $(objdir)run_transaction.o $(openmp)
#objects in scripts/ end


#objects in kvstore/ begin

#objects in sitree/ begin
$(objdir)SITree.o: KVstore/SITree/SITree.cpp KVstore/SITree/SITree.h $(objdir)Stream.o
	$(CXX) $(CFLAGS) KVstore/SITree/SITree.cpp -o $(objdir)SITree.o $(openmp)

$(objdir)SIStorage.o: KVstore/SITree/storage/SIStorage.cpp KVstore/SITree/storage/SIStorage.h $(objdir)Util.o
	$(CXX) $(CFLAGS) KVstore/SITree/storage/SIStorage.cpp -o $(objdir)SIStorage.o $(def64IO) $(openmp)

$(objdir)SINode.o: KVstore/SITree/node/SINode.cpp KVstore/SITree/node/SINode.h $(objdir)Util.o
	$(CXX) $(CFLAGS) KVstore/SITree/node/SINode.cpp -o $(objdir)SINode.o $(openmp)

$(objdir)SIIntlNode.o: KVstore/SITree/node/SIIntlNode.cpp KVstore/SITree/node/SIIntlNode.h
	$(CXX) $(CFLAGS) KVstore/SITree/node/SIIntlNode.cpp -o $(objdir)SIIntlNode.o $(openmp)

$(objdir)SILeafNode.o: KVstore/SITree/node/SILeafNode.cpp KVstore/SITree/node/SILeafNode.h
	$(CXX) $(CFLAGS) KVstore/SITree/node/SILeafNode.cpp -o $(objdir)SILeafNode.o $(openmp)

$(objdir)SIHeap.o: KVstore/SITree/heap/SIHeap.cpp KVstore/SITree/heap/SIHeap.h $(objdir)Util.o
	$(CXX) $(CFLAGS) KVstore/SITree/heap/SIHeap.cpp -o $(objdir)SIHeap.o $(openmp)
#objects in sitree/ end

#objects in istree/ begin
$(objdir)ISTree.o: KVstore/ISTree/ISTree.cpp KVstore/ISTree/ISTree.h $(objdir)Stream.o
	$(CXX) $(CFLAGS) KVstore/ISTree/ISTree.cpp -o $(objdir)ISTree.o $(openmp)

$(objdir)ISStorage.o: KVstore/ISTree/storage/ISStorage.cpp KVstore/ISTree/storage/ISStorage.h $(objdir)Util.o
	$(CXX) $(CFLAGS) KVstore/ISTree/storage/ISStorage.cpp -o $(objdir)ISStorage.o $(def64IO) $(openmp)

$(objdir)ISNode.o: KVstore/ISTree/node/ISNode.cpp KVstore/ISTree/node/ISNode.h $(objdir)Util.o
	$(CXX) $(CFLAGS) KVstore/ISTree/node/ISNode.cpp -o $(objdir)ISNode.o $(openmp)

$(objdir)ISIntlNode.o: KVstore/ISTree/node/ISIntlNode.cpp KVstore/ISTree/node/ISIntlNode.h
	$(CXX) $(CFLAGS) KVstore/ISTree/node/ISIntlNode.cpp -o $(objdir)ISIntlNode.o $(openmp)

$(objdir)ISLeafNode.o: KVstore/ISTree/node/ISLeafNode.cpp KVstore/ISTree/node/ISLeafNode.h
	$(CXX) $(CFLAGS) KVstore/ISTree/node/ISLeafNode.cpp -o $(objdir)ISLeafNode.o $(openmp)

$(objdir)ISHeap.o: KVstore/ISTree/heap/ISHeap.cpp KVstore/ISTree/heap/ISHeap.h $(objdir)Util.o
	$(CXX) $(CFLAGS) KVstore/ISTree/heap/ISHeap.cpp -o $(objdir)ISHeap.o $(openmp)
#objects in istree/ end

#objects in isarray/ begin
$(objdir)ISArray.o: KVstore/ISArray/ISArray.cpp KVstore/ISArray/ISArray.h $(objdir)VList.o
	$(CXX) $(CFLAGS) KVstore/ISArray/ISArray.cpp -o $(objdir)ISArray.o

$(objdir)ISBlockManager.o: KVstore/ISArray/ISBlockManager.cpp KVstore/ISArray/ISBlockManager.h 
	$(CXX) $(CFLAGS) KVstore/ISArray/ISBlockManager.cpp -o $(objdir)ISBlockManager.o

$(objdir)ISEntry.o: KVstore/ISArray/ISEntry.cpp KVstore/ISArray/ISEntry.h
	$(CXX) $(CFLAGS) KVstore/ISArray/ISEntry.cpp -o $(objdir)ISEntry.o
#objects in isarray/ end

#objects in ivtree/ begin
$(objdir)IVTree.o: KVstore/IVTree/IVTree.cpp KVstore/IVTree/IVTree.h $(objdir)Stream.o $(objdir)VList.o
	$(CXX) $(CFLAGS) KVstore/IVTree/IVTree.cpp -o $(objdir)IVTree.o $(openmp)

$(objdir)IVStorage.o: KVstore/IVTree/storage/IVStorage.cpp KVstore/IVTree/storage/IVStorage.h $(objdir)Util.o
	$(CXX) $(CFLAGS) KVstore/IVTree/storage/IVStorage.cpp -o $(objdir)IVStorage.o $(def64IO) $(openmp)

$(objdir)IVNode.o: KVstore/IVTree/node/IVNode.cpp KVstore/IVTree/node/IVNode.h $(objdir)Util.o
	$(CXX) $(CFLAGS) KVstore/IVTree/node/IVNode.cpp -o $(objdir)IVNode.o $(openmp)

$(objdir)IVIntlNode.o: KVstore/IVTree/node/IVIntlNode.cpp KVstore/IVTree/node/IVIntlNode.h
	$(CXX) $(CFLAGS) KVstore/IVTree/node/IVIntlNode.cpp -o $(objdir)IVIntlNode.o $(openmp)

$(objdir)IVLeafNode.o: KVstore/IVTree/node/IVLeafNode.cpp KVstore/IVTree/node/IVLeafNode.h
	$(CXX) $(CFLAGS) KVstore/IVTree/node/IVLeafNode.cpp -o $(objdir)IVLeafNode.o $(openmp)

$(objdir)IVHeap.o: KVstore/IVTree/heap/IVHeap.cpp KVstore/IVTree/heap/IVHeap.h $(objdir)Util.o
	$(CXX) $(CFLAGS) KVstore/IVTree/heap/IVHeap.cpp -o $(objdir)IVHeap.o $(openmp)
#objects in ivtree/ end

#objects in ivarray/ begin
$(objdir)IVArray.o: KVstore/IVArray/IVArray.cpp KVstore/IVArray/IVArray.h $(objdir)VList.o  $(objdir)Transaction.o 
	$(CXX) $(CFLAGS) KVstore/IVArray/IVArray.cpp -o $(objdir)IVArray.o

$(objdir)IVBlockManager.o: KVstore/IVArray/IVBlockManager.cpp KVstore/IVArray/IVBlockManager.h 
	$(CXX) $(CFLAGS) KVstore/IVArray/IVBlockManager.cpp -o $(objdir)IVBlockManager.o

$(objdir)IVEntry.o: KVstore/IVArray/IVEntry.cpp KVstore/IVArray/IVEntry.h $(objdir)Version.o
	$(CXX) $(CFLAGS) KVstore/IVArray/IVEntry.cpp -o $(objdir)IVEntry.o

#objects in ivarray/ end

$(objdir)KVstore.o: KVstore/KVstore.cpp KVstore/KVstore.h KVstore/Tree.h 
	$(CXX) $(CFLAGS) KVstore/KVstore.cpp $(inc) -o $(objdir)KVstore.o $(openmp)

#objects in kvstore/ end


#objects in Database/ begin

$(objdir)Statistics.o: Database/Statistics.cpp Database/Statistics.h $(objdir)Util.o $(objdir)KVstore.o
	$(CXX) $(CFLAGS) Database/Statistics.cpp $(inc) -o $(objdir)Statistics.o $(openmp)

$(objdir)Database.o: Database/Database.cpp Database/Database.h \
	$(objdir)IDList.o $(objdir)ResultSet.o $(objdir)SPARQLquery.o \
	$(objdir)BasicQuery.o $(objdir)Triple.o $(objdir)SigEntry.o \
	$(objdir)KVstore.o $(objdir)VSTree.o $(objdir)DBparser.o $(objdir)Statistics.o $(objdir)HopIndex.o\
	$(objdir)Util.o $(objdir)RDFParser.o $(objdir)Join.o $(objdir)GeneralEvaluation.o $(objdir)StringIndex.o $(objdir)Transaction.o
	$(CXX) $(CFLAGS) Database/Database.cpp $(inc) -o $(objdir)Database.o $(openmp)

$(objdir)Join.o: Database/Join.cpp Database/Join.h $(objdir)IDList.o $(objdir)BasicQuery.o $(objdir)Util.o\
	$(objdir)KVstore.o $(objdir)Util.o $(objdir)SPARQLquery.o $(objdir)Transaction.o
	$(CXX) $(CFLAGS) Database/Join.cpp $(inc) -o $(objdir)Join.o $(openmp)

$(objdir)Strategy.o: Database/Strategy.cpp Database/Strategy.h $(objdir)SPARQLquery.o $(objdir)BasicQuery.o \
	$(objdir)Triple.o $(objdir)IDList.o $(objdir)KVstore.o $(objdir)VSTree.o $(objdir)Util.o $(objdir)Join.o $(objdir)Transaction.o
	$(CXX) $(CFLAGS) Database/Strategy.cpp $(inc) -o $(objdir)Strategy.o $(openmp)

$(objdir)CSR.o: Database/CSR.cpp Database/CSR.h $(objdir)Util.o
	$(CXX) $(CFLAGS) Database/CSR.cpp $(inc) -o $(objdir)CSR.o $(openmp)

$(objdir)TableOperator.o: Database/TableOperator.cpp Database/TableOperator.h $(objdir)Util.o $(objdir)BasicQuery.o
	$(CXX) $(CFLAGS) Database/TableOperator.cpp $(inc) -o $(objdir)TableOperator.o $(openmp)

$(objdir)ResultTrigger.o: Database/ResultTrigger.cpp Database/ResultTrigger.h $(objdir)Util.o
	$(CXX) $(CFLAGS) Database/ResultTrigger.cpp $(inc) -o $(objdir)ResultTrigger.o $(openmp)

$(objdir)PlanTree.o: Database/PlanTree.cpp Database/PlanTree.h $(objdir)BasicQuery.o $(objdir)TableOperator.o \
    $(objdir)BGPQuery.o $(objdir)Util.o
	$(CXX) $(CFLAGS) Database/PlanTree.cpp $(inc) -o $(objdir)PlanTree.o  $(openmp)

$(objdir)PlanGenerator.o: Database/PlanGenerator.cpp Database/PlanGenerator.h \
	$(objdir)Util.o $(objdir)BasicQuery.o $(objdir)BGPQuery.o $(objdir)IDList.o $(objdir)KVstore.o \
	$(objdir)Statistics.o $(objdir)PlanTree.o $(objdir)TableOperator.o $(objdir)OrderedVector.o
	$(CXX) $(CFLAGS) Database/PlanGenerator.cpp $(inc) -o $(objdir)PlanGenerator.o $(openmp)

$(objdir)Executor.o: Database/Executor.cpp Database/Executor.h $(objdir)Util.o $(objdir)SPARQLquery.o $(objdir)BasicQuery.o $(objdir)IDList.o \
	$(objdir)KVstore.o  $(objdir)VSTree.o $(objdir)Join.o $(objdir)Transaction.o $(objdir)TableOperator.o $(objdir)ResultTrigger.o \
	$(objdir)QueryPlan.o $(objdir)Statistics.o $(objdir)PlanTree.o $(objdir)PlanGenerator.o $(objdir)OrderedVector.o \
	$(objdir)DPBTopKUtil.o $(objdir)DPPTopKUtil.o $(objdir)kTPMUtil.o $(objdir)SAEUtil.o $(objdir)DfsUtil.o \
	$(objdir)DfsUtilOld.o $(objdir)EagerUtil.o $(objdir)Take2Util.o $(objdir)TakeAllUtil.o
	$(CXX) $(CFLAGS) Database/Executor.cpp $(inc) -o $(objdir)Executor.o $(openmp)

$(objdir)Optimizer.o: Database/Optimizer.cpp Database/Optimizer.h $(objdir)Util.o $(objdir)SPARQLquery.o $(objdir)BasicQuery.o $(objdir)IDList.o \
	$(objdir)KVstore.o  $(objdir)VSTree.o $(objdir)Join.o $(objdir)Transaction.o $(objdir)TableOperator.o $(objdir)ResultTrigger.o $(objdir)Executor.o\
	$(objdir)QueryPlan.o $(objdir)Statistics.o $(objdir)PlanTree.o $(objdir)PlanGenerator.o $(objdir)OrderedVector.o \
	$(objdir)DPBTopKUtil.o $(objdir)DPPTopKUtil.o $(objdir)kTPMUtil.o $(objdir)SAEUtil.o $(objdir)DfsUtil.o \
	$(objdir)DfsUtilOld.o $(objdir)EagerUtil.o $(objdir)Take2Util.o $(objdir)TakeAllUtil.o
	$(CXX) $(CFLAGS) Database/Optimizer.cpp $(inc) -o $(objdir)Optimizer.o $(openmp)

$(objdir)Txn_manager.o: Database/Txn_manager.cpp Database/Txn_manager.h $(objdir)Util.o $(objdir)Transaction.o $(objdir)Database.o
	$(CXX) $(CFLAGS) Database/Txn_manager.cpp $(inc) -o $(objdir)Txn_manager.o $(openmp)

#objects in Database/ end


#objects in Query/ begin

$(objdir)IDList.o: Query/IDList.cpp Query/IDList.h
	$(CXX) $(CFLAGS) Query/IDList.cpp $(inc) -o $(objdir)IDList.o $(openmp)

$(objdir)SPARQLquery.o: Query/SPARQLquery.cpp Query/SPARQLquery.h $(objdir)BasicQuery.o
	$(CXX) $(CFLAGS) Query/SPARQLquery.cpp $(inc) -o $(objdir)SPARQLquery.o $(openmp)

$(objdir)BasicQuery.o: Query/BasicQuery.cpp Query/BasicQuery.h $(objdir)Signature.o
	$(CXX) $(CFLAGS) Query/BasicQuery.cpp $(inc) -o $(objdir)BasicQuery.o $(openmp)

$(objdir)ResultSet.o: Query/ResultSet.cpp Query/ResultSet.h $(objdir)Stream.o
	$(CXX) $(CFLAGS) Query/ResultSet.cpp $(inc) -o $(objdir)ResultSet.o $(openmp)

$(objdir)Varset.o: Query/Varset.cpp Query/Varset.h
	$(CXX) $(CFLAGS) Query/Varset.cpp $(inc) -o $(objdir)Varset.o $(openmp)

$(objdir)QueryPlan.o: Query/QueryPlan.cpp Query/QueryPlan.h $(objdir)BasicQuery.o $(objdir)TableOperator.o $(objdir)ResultTrigger.o
	$(CXX) $(CFLAGS) Query/QueryPlan.cpp $(inc) -o $(objdir)QueryPlan.o $(openmp)

$(objdir)QueryTree.o: Query/QueryTree.cpp Query/QueryTree.h $(objdir)Varset.o
	$(CXX) $(CFLAGS) Query/QueryTree.cpp $(inc) -o $(objdir)QueryTree.o $(openmp)

$(objdir)TempResult.o: Query/TempResult.cpp Query/TempResult.h Query/RegexExpression.h $(objdir)Util.o \
	$(objdir)StringIndex.o $(objdir)QueryTree.o $(objdir)Varset.o $(objdir)EvalMultitypeValue.o
	$(CXX) $(CFLAGS) Query/TempResult.cpp $(inc) -o $(objdir)TempResult.o $(openmp)

$(objdir)QueryCache.o: Query/QueryCache.cpp Query/QueryCache.h $(objdir)Util.o $(objdir)QueryTree.o \
	$(objdir)TempResult.o $(objdir)Varset.o
	$(CXX) $(CFLAGS) Query/QueryCache.cpp $(inc) -o $(objdir)QueryCache.o $(openmp)

$(objdir)PathQueryHandler.o: Query/PathQueryHandler.cpp Query/PathQueryHandler.h $(objdir)Util.o $(objdir)CSR.o
	$(CXX) $(CFLAGS) Query/PathQueryHandler.cpp $(inc) -o $(objdir)PathQueryHandler.o $(openmp)

$(objdir)BGPQuery.o: Query/BGPQuery.cpp Query/BGPQuery.h $(objdir)BasicQuery.o  $(objdir)Util.o \
    $(objdir)Triple.o $(objdir)KVstore.o
	$(CXX) $(CFLAGS) Query/BGPQuery.cpp $(inc) -o $(objdir)BGPQuery.o $(openmp)

#objects in Query/topk/ begin

$(objdir)HopIndex.o: Query/topk/PTAB Query/topk/PTAB $(objdir)Util.o $(objdir)KVstore.o
	$(CXX) $(CFLAGS) Query/topk/PTAB/HopIndex.cpp $(inc) -o $(objdir)HopIndex.o $(openmp)

$(objdir)Pool.o: Query/topk/CommonDP/Pool.cpp Query/topk/CommonDP/Pool.h $(objdir)Util.o
	$(CXX) $(CFLAGS) Query/topk/CommonDP/Pool.cpp $(inc) -o $(objdir)Pool.o $(openmp)

$(objdir)DynamicTrie.o: Query/topk/CommonDP/DynamicTrie.cpp Query/topk/CommonDP/DynamicTrie.h $(objdir)Util.o $(objdir)Pool.o
	$(CXX) $(CFLAGS) Query/topk/CommonDP/DynamicTrie.cpp $(inc) -o $(objdir)DynamicTrie.o $(openmp)

$(objdir)DPBList.o: Query/topk/DPB/DPBList.cpp Query/topk/DPB/DPBList.h $(objdir)Util.o $(objdir)Pool.o \
	$(objdir)DynamicTrie.o
	$(CXX) $(CFLAGS) Query/topk/DPB/DPBList.cpp $(inc) -o $(objdir)DPBList.o $(openmp)

$(objdir)DPPList.o: Query/topk/DPP/DPPList.cpp Query/topk/DPP/DPPList.h $(objdir)Util.o $(objdir)Pool.o \
	$(objdir)DynamicTrie.o
	$(CXX) $(CFLAGS) Query/topk/DPP/DPPList.cpp $(inc) -o $(objdir)DPPList.o $(openmp)

$(objdir)GlobalQueue.o: Query/topk/DPP/GlobalQueue.cpp Query/topk/DPP/GlobalQueue.h $(objdir)DPPList.o
	$(CXX) $(CFLAGS) Query/topk/DPP/GlobalQueue.cpp $(inc) -o $(objdir)GlobalQueue.o $(openmp)

$(objdir)DPPFQIterator.o: Query/topk/DPP/DPPFQIterator.cpp Query/topk/DPP/DPPFQIterator.h $(objdir)Util.o $(objdir)Pool.o \
	$(objdir)DynamicTrie.o $(objdir)DPPFRIterator.o $(objdir)DPPOWIterator.o $(objdir)GlobalQueue.o
	$(CXX) $(CFLAGS) Query/topk/DPP/DPPFQIterator.cpp $(inc) -o $(objdir)DPPFQIterator.o $(openmp)

$(objdir)DPPFRIterator.o: Query/topk/DPP/DPPFRIterator.cpp Query/topk/DPP/DPPFRIterator.h $(objdir)Util.o $(objdir)Pool.o \
	$(objdir)DynamicTrie.o $(objdir)DPPFQIterator.o $(objdir)DPPOWIterator.o $(objdir)GlobalQueue.o
	$(CXX) $(CFLAGS) Query/topk/DPP/DPPFRIterator.cpp $(inc) -o $(objdir)DPPFRIterator.o $(openmp)

$(objdir)DPPOWIterator.o: Query/topk/DPP/DPPOWIterator.cpp Query/topk/DPP/DPPOWIterator.h $(objdir)Util.o $(objdir)Pool.o \
	$(objdir)DynamicTrie.o $(objdir)DPPFRIterator.o $(objdir)DPPFQIterator.o $(objdir)GlobalQueue.o
	$(CXX) $(CFLAGS) Query/topk/DPP/DPPOWIterator.cpp $(inc) -o $(objdir)DPPOWIterator.o $(openmp)

$(objdir)Subspace.o: Query/topk/kTPM/Subspace.cpp Query/topk/kTPM/Subspace.h $(objdir)Util.o $(objdir)Pool.o
	$(CXX) $(CFLAGS) Query/topk/kTPM/Subspace.cpp $(inc) -o $(objdir)Subspace.o $(openmp)

$(objdir)SpaceSAE.o: Query/topk/SAE/SpaceSAE.cpp Query/topk/SAE/SpaceSAE.h $(objdir)Util.o $(objdir)Pool.o
	$(CXX) $(CFLAGS) Query/topk/SAE/SpaceSAE.cpp $(inc) -o $(objdir)SpaceSAE.o $(openmp)

$(objdir)EagerSubspace.o: Query/topk/EAGER/EagerSubspace.cpp Query/topk/EAGER/EagerSubspace.h $(objdir)Util.o $(objdir)Pool.o
	$(CXX) $(CFLAGS) Query/topk/EAGER/EagerSubspace.cpp $(inc) -o $(objdir)EagerSubspace.o $(openmp)

$(objdir)Take2Subspace.o: Query/topk/TAKE2/Take2Subspace.cpp Query/topk/TAKE2/Take2Subspace.h $(objdir)Util.o $(objdir)Pool.o
	$(CXX) $(CFLAGS) Query/topk/TAKE2/Take2Subspace.cpp $(inc) -o $(objdir)Take2Subspace.o $(openmp)

$(objdir)TakeAllSubspace.o: Query/topk/TAKEALL/TakeAllSubspace.cpp Query/topk/TAKEALL/TakeAllSubspace.h $(objdir)Util.o $(objdir)Pool.o
	$(CXX) $(CFLAGS) Query/topk/TAKEALL/TakeAllSubspace.cpp $(inc) -o $(objdir)TakeAllSubspace.o $(openmp)


$(objdir)TopKSearchPlan.o: Query/topk/TopKSearchPlan.cpp Query/topk/TopKSearchPlan.h $(objdir)Util.o \
	$(objdir)KVstore.o $(objdir)SPARQLquery.o $(objdir)BasicQuery.o \
	$(objdir)Statistics.o $(objdir)QueryTree.o $(objdir)IDList.o \
	$(objdir)PlanGenerator.o $(objdir)TableOperator.o
	$(CXX) $(CFLAGS) Query/topk/TopKSearchPlan.cpp $(inc) -o $(objdir)TopKSearchPlan.o $(openmp)

$(objdir)TopKUtil.o: Query/topk/TopKUtil.cpp Query/topk/TopKUtil.h $(objdir)Util.o $(objdir)KVstore.o \
	$(objdir)SPARQLquery.o $(objdir)BasicQuery.o $(objdir)Statistics.o \
	$(objdir)QueryTree.o $(objdir)IDList.o $(objdir)TableOperator.o $(objdir)TopKSearchPlan.o
	$(CXX) $(CFLAGS) Query/topk/TopKUtil.cpp $(inc) -o $(objdir)TopKUtil.o $(openmp)

$(objdir)DPBTopKUtil.o: Query/topk/DPB/DPBTopKUtil.cpp Query/topk/DPB/DPBTopKUtil.h $(objdir)TopKUtil.o $(objdir)DPBList.o
	$(CXX) $(CFLAGS) Query/topk/DPB/DPBTopKUtil.cpp $(inc) -o $(objdir)DPBTopKUtil.o $(openmp)

$(objdir)DPPTopKUtil.o: Query/topk/DPP/DPPTopKUtil.cpp Query/topk/DPP/DPPTopKUtil.h $(objdir)TopKUtil.o \
	$(objdir)DPPFQIterator.o  $(objdir)DPPFRIterator.o  $(objdir)DPPOWIterator.o
	$(CXX) $(CFLAGS) Query/topk/DPP/DPPTopKUtil.cpp $(inc) -o $(objdir)DPPTopKUtil.o $(openmp)

$(objdir)kTPMUtil.o: Query/topk/kTPM/kTPMUtil.cpp Query/topk/kTPM/kTPMUtil.h $(objdir)TopKUtil.o $(objdir)Subspace.o
	$(CXX) $(CFLAGS) Query/topk/kTPM/kTPMUtil.cpp $(inc) -o $(objdir)kTPMUtil.o $(openmp)

$(objdir)SAEUtil.o: Query/topk/SAE/SAEUtil.cpp Query/topk/SAE/SAEUtil.h $(objdir)TopKUtil.o $(objdir)Subspace.o
	$(CXX) $(CFLAGS) Query/topk/SAE/SAEUtil.cpp  $(inc) -o $(objdir)SAEUtil.o $(openmp)

$(objdir)EagerUtil.o: Query/topk/EAGER/EagerUtil.cpp Query/topk/EAGER/EagerUtil.h $(objdir)TopKUtil.o $(objdir)EagerSubspace.o
	$(CXX) $(CFLAGS) Query/topk/EAGER/EagerUtil.cpp $(inc) -o $(objdir)EagerUtil.o $(openmp)

$(objdir)Take2Util.o: Query/topk/TAKE2/Take2Util.cpp Query/topk/TAKE2/Take2Util.h $(objdir)TopKUtil.o $(objdir)Take2Subspace.o
	$(CXX) $(CFLAGS) Query/topk/TAKE2/Take2Util.cpp $(inc) -o $(objdir)Take2Util.o $(openmp)

$(objdir)TakeAllUtil.o: Query/topk/TAKEALL/TakeAllUtil.cpp Query/topk/TAKEALL/TakeAllUtil.h $(objdir)TopKUtil.o $(objdir)TakeAllSubspace.o
	$(CXX) $(CFLAGS) Query/topk/TAKEALL/TakeAllUtil.cpp $(inc) -o $(objdir)TakeAllUtil.o $(openmp)

$(objdir)CompressedVector.o: Query/topk/PTAB Query/topk/PTAB $(objdir)Util.o
	$(CXX) $(CFLAGS) Query/topk/PTAB/CompressedVector.cpp $(inc) -o $(objdir)CompressedVector.o $(openmp)

$(objdir)DfsList.o: Query/topk/PTAB Query/topk/PTAB $(objdir)Util.o $(objdir)Pool.o \
	$(objdir)CompressedVector.o $(objdir)IDList.o $(objdir)TopKSearchPlan.o $(objdir)TopKUtil.o
	$(CXX) $(CFLAGS) Query/topk/PTAB/DfsList.cpp $(inc) -o $(objdir)DfsList.o $(openmp)

$(objdir)DfsFQIterator.o: Query/topk/PTAB Query/topk/PTAB $(objdir)Util.o $(objdir)Pool.o \
	$(objdir)CompressedVector.o $(objdir)DfsFRIterator.o $(objdir)DfsOWIterator.o
	$(CXX) $(CFLAGS) Query/topk/PTAB/DfsFQIterator.cpp $(inc) -o $(objdir)DfsFQIterator.o $(openmp)

$(objdir)DfsFQIteratorOld.o: Query/topk/PTAB Query/topk/PTAB $(objdir)Util.o $(objdir)Pool.o \
	$(objdir)DfsFRIterator.o $(objdir)DfsOWIterator.o $(objdir)DynamicTrie.o
	$(CXX) $(CFLAGS) Query/topk/PTAB/DfsFQIteratorOld.cpp $(inc) -o $(objdir)DfsFQIteratorOld.o $(openmp)

$(objdir)DfsFRIterator.o:  Query/topk/PTAB Query/topk/PTAB $(objdir)Util.o $(objdir)Pool.o \
	$(objdir)CompressedVector.o $(objdir)DfsFQIterator.o $(objdir)DfsOWIterator.o $(objdir)DfsFQIteratorOld.o
	$(CXX) $(CFLAGS) Query/topk/PTAB/DfsFRIterator.cpp $(inc) -o $(objdir)DfsFRIterator.o $(openmp)

$(objdir)DfsOWIterator.o: Query/topk/PTAB Query/topk/PTAB $(objdir)Util.o $(objdir)Pool.o \
	$(objdir)CompressedVector.o $(objdir)DfsFRIterator.o $(objdir)DfsFQIterator.o $(objdir)DfsFQIteratorOld.o
	$(CXX) $(CFLAGS) Query/topk/PTAB/DfsOWIterator.cpp $(inc) -o $(objdir)DfsOWIterator.o $(openmp)

$(objdir)DfsUtil.o: Query/topk/PTAB Query/topk/PTAB $(objdir)TopKUtil.o \
	$(objdir)DfsFQIterator.o  $(objdir)DfsFRIterator.o  $(objdir)DfsOWIterator.o  $(objdir)Pool.o
	$(CXX) $(CFLAGS) Query/topk/PTAB/DfsUtil.cpp $(inc) -o $(objdir)DfsUtil.o $(openmp)

$(objdir)DfsUtilOld.o: Query/topk/PTAB Query/topk/PTAB $(objdir)TopKUtil.o \
	$(objdir)DfsFQIteratorOld.o  $(objdir)DfsFRIterator.o  $(objdir)DfsOWIterator.o  $(objdir)Pool.o
	$(CXX) $(CFLAGS) Query/topk/PTAB/DfsUtilOld.cpp $(inc) -o $(objdir)DfsUtilOld.o $(openmp)
#objects in Query/topk/ end


#no more using $(objdir)Database.o
$(objdir)GeneralEvaluation.o: Query/GeneralEvaluation.cpp Query/GeneralEvaluation.h Query/RegexExpression.h \
	$(objdir)VSTree.o $(objdir)KVstore.o $(objdir)StringIndex.o $(objdir)Strategy.o $(objdir)QueryParser.o \
	$(objdir)Triple.o $(objdir)Util.o $(objdir)EvalMultitypeValue.o $(objdir)SPARQLquery.o $(objdir)QueryTree.o $(objdir)Varset.o  $(objdir)Statistics.o \
	$(objdir)TempResult.o $(objdir)QueryCache.o $(objdir)ResultSet.o $(objdir)PathQueryHandler.o $(objdir)Optimizer.o
	$(CXX) $(CFLAGS) Query/GeneralEvaluation.cpp $(inc) -o $(objdir)GeneralEvaluation.o $(openmp)

#objects in Query/ end


#objects in Signature/ begin

$(objdir)SigEntry.o: Signature/SigEntry.cpp Signature/SigEntry.h $(objdir)Signature.o
	$(CXX) $(CFLAGS) Signature/SigEntry.cpp $(inc) -o $(objdir)SigEntry.o $(openmp)

$(objdir)Signature.o: Signature/Signature.cpp Signature/Signature.h
	$(CXX) $(CFLAGS) Signature/Signature.cpp $(inc) -o $(objdir)Signature.o $(openmp)

#objects in Signature/ end


#objects in Util/ begin

$(objdir)Util.o:  Util/Util.cpp Util/Util.h
	$(CXX) $(CFLAGS) Util/Util.cpp -o $(objdir)Util.o $(openmp)

$(objdir)Stream.o:  Util/Stream.cpp Util/Stream.h $(objdir)Util.o $(objdir)Bstr.o
	$(CXX) $(CFLAGS) Util/Stream.cpp -o $(objdir)Stream.o $(def64IO) $(openmp)

$(objdir)Bstr.o: Util/Bstr.cpp Util/Bstr.h $(objdir)Util.o
	$(CXX) $(CFLAGS)  Util/Bstr.cpp -o $(objdir)Bstr.o $(openmp)

$(objdir)Triple.o: Util/Triple.cpp Util/Triple.h $(objdir)Util.o
	$(CXX) $(CFLAGS) Util/Triple.cpp -o $(objdir)Triple.o $(openmp)

$(objdir)BloomFilter.o:  Util/BloomFilter.cpp Util/BloomFilter.h $(objdir)Util.o
	$(CXX) $(CFLAGS) Util/BloomFilter.cpp -o $(objdir)BloomFilter.o $(openmp) 

$(objdir)VList.o:  Util/VList.cpp Util/VList.h
	$(CXX) $(CFLAGS) Util/VList.cpp -o $(objdir)VList.o $(openmp)

$(objdir)EvalMultitypeValue.o: Util/EvalMultitypeValue.cpp Util/EvalMultitypeValue.h
	$(CXX) $(CFLAGS) Util/EvalMultitypeValue.cpp -o $(objdir)EvalMultitypeValue.o $(openmp)

$(objdir)Version.o: Util/Version.cpp Util/Version.h
	$(CXX) $(CFLAGS) Util/Version.cpp -o $(objdir)Version.o $(openmp)

$(objdir)Transaction.o: Util/Transaction.cpp Util/Transaction.h $(objdir)Util.o $(objdir)IDTriple.o
	$(CXX) $(CFLAGS) Util/Transaction.cpp $(inc) -o $(objdir)Transaction.o $(openmp)

$(objdir)IDTriple.o: Util/IDTriple.cpp Util/IDTriple.h
	$(CXX) $(CFLAGS) Util/IDTriple.cpp -o $(objdir)IDTriple.o $(openmp)

$(objdir)Latch.o: Util/Latch.cpp Util/Latch.h
	$(CXX) $(CFLAGS) Util/Latch.cpp -o $(objdir)Latch.o $(openmp)

$(objdir)IPWhiteList.o:  Util/IPWhiteList.cpp Util/IPWhiteList.h $(objdir)Util.o
	$(CXX) $(CFLAGS) Util/IPWhiteList.cpp -o $(objdir)IPWhiteList.o $(def64IO) $(openmp)

$(objdir)IPBlackList.o:  Util/IPBlackList.cpp Util/IPBlackList.h $(objdir)Util.o
	$(CXX) $(CFLAGS) Util/IPBlackList.cpp -o $(objdir)IPBlackList.o $(def64IO) $(openmp)

$(objdir)OrderedVector.o: Util/OrderedVector.cpp Util/OrderedVector.h
	$(CXX) $(CFLAGS) Util/OrderedVector.cpp -o $(objdir)OrderedVector.o $(openmp)

#objects in util/ end


#objects in VSTree/ begin

$(objdir)VSTree.o: VSTree/VSTree.cpp VSTree/VSTree.h $(objdir)EntryBuffer.o $(objdir)LRUCache.o $(objdir)VNode.o
	$(CXX) $(CFLAGS) VSTree/VSTree.cpp $(inc) -o $(objdir)VSTree.o $(def64IO) $(openmp)

$(objdir)EntryBuffer.o: VSTree/EntryBuffer.cpp VSTree/EntryBuffer.h Signature/SigEntry.h
	$(CXX) $(CFLAGS) VSTree/EntryBuffer.cpp $(inc) -o $(objdir)EntryBuffer.o $(def64IO) $(openmp)

$(objdir)LRUCache.o: VSTree/LRUCache.cpp  VSTree/LRUCache.h VSTree/VNode.h
	$(CXX) $(CFLAGS) VSTree/LRUCache.cpp $(inc) -o $(objdir)LRUCache.o $(def64IO) $(openmp)

$(objdir)VNode.o: VSTree/VNode.cpp VSTree/VNode.h
	$(CXX) $(CFLAGS) VSTree/VNode.cpp $(inc) -o $(objdir)VNode.o $(def64IO) $(openmp)

#objects in VSTree/ end


#objects in StringIndex/ begin
$(objdir)StringIndex.o: StringIndex/StringIndex.cpp StringIndex/StringIndex.h $(objdir)KVstore.o $(objdir)Util.o
	$(CXX) $(CFLAGS) StringIndex/StringIndex.cpp $(inc) -o $(objdir)StringIndex.o $(def64IO) $(openmp)
#objects in StringIndex/ end


#objects in Parser/ begin

$(objdir)DBparser.o: Parser/DBparser.cpp Parser/DBparser.h $(objdir)SPARQLParser.o $(objdir)SPARQLLexer.o $(objdir)Triple.o
	$(CXX) $(CFLAGS) Parser/DBparser.cpp $(inc) -o $(objdir)DBparser.o $(openmp)

# $(objdir)SparqlParser.o: Parser/SparqlParser.c Parser/SparqlParser.h
# 	gcc  $(CFLAGS) Parser/SparqlParser.c $(inc) -o $(objdir)SparqlParser.o $(openmp)

$(objdir)SPARQLParser.o: Parser/SPARQL/SPARQLParser.cpp Parser/SPARQL/SPARQLParser.h
	$(CXX)  $(CFLAGS) Parser/SPARQL/SPARQLParser.cpp $(inc) -o $(objdir)SPARQLParser.o $(openmp)

# $(objdir)SparqlLexer.o: Parser/SparqlLexer.c Parser/SparqlLexer.h
# 	gcc  $(CFLAGS) Parser/SparqlLexer.c $(inc) -o $(objdir)SparqlLexer.o $(openmp)

$(objdir)SPARQLLexer.o: Parser/SPARQL/SPARQLLexer.cpp Parser/SPARQL/SPARQLLexer.h
	$(CXX)  $(CFLAGS) Parser/SPARQL/SPARQLLexer.cpp $(inc) -o $(objdir)SPARQLLexer.o $(openmp)

$(objdir)TurtleParser.o: Parser/TurtleParser.cpp Parser/TurtleParser.h Parser/Type.h
	gcc  $(CFLAGS) Parser/TurtleParser.cpp $(inc) -o $(objdir)TurtleParser.o $(openmp)

$(objdir)RDFParser.o: Parser/RDFParser.cpp Parser/RDFParser.h $(objdir)TurtleParser.o $(objdir)Triple.o
	gcc  $(CFLAGS) Parser/RDFParser.cpp $(inc) -o $(objdir)RDFParser.o $(openmp)

$(objdir)QueryParser.o: Parser/QueryParser.cpp Parser/QueryParser.h $(objdir)SPARQLParser.o $(objdir)SPARQLLexer.o $(objdir)QueryTree.o
	$(CXX) $(CFLAGS) Parser/QueryParser.cpp $(inc) -o $(objdir)QueryParser.o $(openmp)

#objects in Parser/ end

#objects in Trie/ begin

$(objdir)TrieNode.o: Trie/TrieNode.cpp Trie/TrieNode.h
	$(CXX) $(CFLAGS) Trie/TrieNode.cpp -o $(objdir)TrieNode.o

$(objdir)Trie.o: Trie/Trie.cpp Trie/Trie.h $(objdir)TrieNode.o $(objdir)Triple.o $(objdir)RDFParser.o
	$(CXX) $(CFLAGS) Trie/Trie.cpp $(inc) -o $(objdir)Trie.o

#objects in Server/ begin

$(objdir)Operation.o: Server/Operation.cpp Server/Operation.h
	$(CXX) $(CFLAGS) Server/Operation.cpp $(inc) -o $(objdir)Operation.o $(openmp)

$(objdir)Socket.o: Server/Socket.cpp Server/Socket.h
	$(CXX) $(CFLAGS) Server/Socket.cpp $(inc) -o $(objdir)Socket.o $(openmp)

$(objdir)Server.o: Server/Server.cpp Server/Server.h $(objdir)Socket.o $(objdir)Database.o $(objdir)Operation.o
	$(CXX) $(CFLAGS) Server/Server.cpp $(inc) -o $(objdir)Server.o $(openmp)

$(objdir)Client.o: Server/Client.cpp Server/Client.h $(objdir)Socket.o $(objdir)Util.o
	$(CXX) $(CFLAGS) Server/Client.cpp $(inc) -o $(objdir)Client.o $(openmp)

# $(objdir)client_http.o: Server/client_http.hpp
# 	$(CXX) $(CFLAGS) Server/client_http.hpp $(inc) -o $(objdir)client_http.o

# $(objdir)server_http.o: Server/server_http.hpp
# 	$(CXX) $(CFLAGS) Server/server_http.hpp $(inc) -o $(objdir)server_http.o

#objects in Server/ end


pre:
	rm -rf tools/rapidjson/
	rm -rf tools/antlr4-cpp-runtime-4
	cd tools; tar -xzvf rapidjson.tar.gz;
	cd tools; tar -xzvf antlr4-cpp-runtime-4.tar.gz;
	cd tools/antlr4-cpp-runtime-4/; cmake .; make; cp dist/libantlr4-runtime.a ../../lib/;

# DEBUG: below not works properly
#Parser/SparqlLexer.c Parser/SparqlLexer.h Parser/SparqlParser.h Parser/SparqlParser.c: unpack_sparql
#.INTERMEDIATE: unpack_sparql
#unpack_sparql: tools/sparql.tar.gz
##NOTICE: update the sparql.tar.gz if Sparql* in Parser are changed manually
#rm -rf Parser/Sparql*
#cd tools; tar -xzvf sparql.tar.gz; mv Sparql* ../Parser/;

$(api_cpp): $(objdir)Socket.o
	$(MAKE) -C api/socket/cpp/src 
	$(MAKE) -C api/http/cpp/

$(api_java):
	$(MAKE) -C api/socket/java/src
	$(MAKE) -C api/http/java/src

.PHONY: clean dist tarball api_example gtest sumlines contribution test

test: $(TARGET)
	@echo "basic build/query/add/sub/drop test......"
	@bash scripts/basic_test.sh
	@echo "repeatedly insertion/deletion test......"
	@scripts/update_test > /dev/null
	@echo "parser test......"
	@bash scripts/parser_test.sh

clean:
	#rm -rf lib/libantlr4-runtime.a
	$(MAKE) -C api/socket/cpp/src clean
	$(MAKE) -C api/socket/cpp/example clean
	$(MAKE) -C api/socket/java/src clean
	$(MAKE) -C api/socket/java/example clean
	$(MAKE) -C api/http/cpp clean
	$(MAKE) -C api/http/cpp/src clean
	$(MAKE) -C api/http/cpp/example clean
	$(MAKE) -C api/http/java/src clean
	$(MAKE) -C api/http/java/example clean
	#$(MAKE) -C KVstore clean
	rm -rf $(exedir)g* $(objdir)*.o $(exedir).gserver* $(exedir)shutdown $(exedir).gconsole* $(exedir)rollback
	rm -rf bin/*.class
	rm -rf $(testdir)update_test $(testdir)dataset_test $(testdir)transaction_test $(testdir)run_transaction
	#rm -rf .project .cproject .settings   just for eclipse
	rm -rf logs/*.log
	rm -rf *.out   # gmon.out for gprof with -pg


dist: clean
	rm -rf *.nt *.n3 .debug/*.log .tmp/*.dat *.txt *.db
	rm -rf lib/libantlr4-runtime.a
	rm -rf cscope* .cproject .settings tags
	rm -rf *.info
	rm -rf backups/*.db

tarball:
	tar -czvf gstore.tar.gz api backups bin lib tools .debug .tmp .objs scripts garbage docs data logs \
		Main Database KVstore Util Query Signature VSTree Parser Server README.md init.conf NOTES.md StringIndex COVERAGE \
		Dockerfile LICENSE makefile Trie

APIexample: $(api_cpp) $(api_java)
	$(MAKE) -C api/socket/cpp/example
	$(MAKE) -C api/socket/java/example
	$(MAKE) -C api/http/cpp/example
	$(MAKE) -C api/http/java/example

gtest: $(objdir)gtest.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)gtest $(objdir)gtest.o $(objfile) lib/libantlr4-runtime.a $(library) $(openmp)

$(objdir)gtest.o: scripts/gtest.cpp
	$(CXX) $(CFLAGS) scripts/gtest.cpp $(inc) -o $(objdir)gtest.o $(openmp)
	
$(exedir)gadd: $(objdir)gadd.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)gadd $(objdir)gadd.o $(objfile) lib/libantlr4-runtime.a $(library) $(openmp)

$(objdir)gadd.o: Main/gadd.cpp
	$(CXX) $(CFLAGS) Main/gadd.cpp $(inc) -o $(objdir)gadd.o $(openmp)

#$(objdir)HttpConnector: $(objdir)HttpConnector.o $(objfile)
	#$(CXX) $(CFLAGS) -o $(exedir)HttpConnector $(objdir)HttpConnector.o $(objfile) lib/libantlr4-runtime.a $(library) $(inc) -DUSE_BOOST_REGEX

#$(objdir)HttpConnector.o: Main/HttpConnector.cpp
	#$(CXX) $(CFLAGS) Main/HttpConnector.cpp $(inc) -o $(objdir)HttpConnector.o $(library) -DUSE_BOOST_REGEX

$(exedir)gsub: $(objdir)gsub.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)gsub $(objdir)gsub.o $(objfile) lib/libantlr4-runtime.a $(library) $(openmp)

$(objdir)gsub.o: Main/gsub.cpp
	$(CXX) $(CFLAGS) Main/gsub.cpp $(inc) -o $(objdir)gsub.o $(openmp)

$(exedir)indexBuilder: $(objdir)indexBuilder.o $(objfile)
	$(CXX) $(EXEFLAG) -o $(exedir)indexBuilder $(objdir)indexBuilder.o $(objfile) lib/libantlr4-runtime.a $(library) $(openmp)

$(objdir)indexBuilder.o: Main/indexBuilder.cpp
	$(CXX) $(CFLAGS) Main/indexBuilder.cpp $(inc) -o $(objdir)indexBuilder.o $(openmp)

sumlines:
	@bash scripts/sumline.sh

tag:
	ctags -R

idx:
	find `realpath .` -name "*.h" -o -name "*.c" -o -name "*.cpp" > cscope.files
	cscope -bkq #-i cscope.files

cover:
	bash scripts/cover.sh

fulltest:
	#NOTICE:compile gstore with -O2 only
	#setup new virtuoso and configure it
	cp scripts/full_test.sh ~
	cd ~
	bash full_test.sh

#test the efficience of kvstore, insert/delete/search, use dbpedia170M by default
test-kvstore:
	# test/kvstore_test.cpp
	echo "TODO"

# https://segmentfault.com/a/1190000008542123
contribution:
	bash scripts/contribution.sh

