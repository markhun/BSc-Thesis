# Specify the path to regent.py
REGENT = /usr/local/legion/language/regent.py #~/Documents/BA_Legion/legion/language/regent.py
# specify the path to legion_prof.py
LEGION_PROF = /usr/local/legion/tools/legion_prof.py #~/Documents/BA_Legion/legion/tools/legion_prof.py
# specify the path to legion_spy.py
LEGION_SPY = /usr/local/legion/tools/legion_spy.py #~/Documents/BA_Legion/legion/tools/legion_spy.py
# specify number of proc to run on:
CPUS = 6

DATE = $(shell date -u +'%Y-%m-%dT%H:%M:%S')
# checking the version of python currently active
# Legion Prof and Spy produce only cryptic errors when you try to run them under Python3
python_version_full := $(wordlist 2,4,$(subst ., ,$(shell python --version 2>&1)))
python_version_major := $(word 1,${python_version_full})

# There are different ways to run a regent program
# just running it
run:
	$(REGENT) sph_partitioning.rg -ll:cpu $(CPUS)
# running it with the Legion Profiler generating data
profiling_data.gz:
	$(REGENT) sph_partitioning.rg -ll:cpu $(CPUS) -lg:prof 1 -lg:prof_logfile profiling_data.gz
# running it with Legion Spy generating data
spy.log:
	$(REGENT) sph_partitioning.rg -ll:cpu $(CPUS) -lg:spy -logfile spy.log
# running it with both the Profiler and Spy - Be careful, this can lead to unexpected side effects!
runspyprof:
	$(REGENT) sph_partitioning.rg -ll:cpu $(CPUS) -lg:prof 1 -lg:prof_logfile profiling_data.gz -lg:spy -logfile spy.log


all: runspyprof
ifeq ($(python_version_major), 3)
	@echo Python 3 Enviroment active. Legion Profiler only supports Python 2.7 !!!
else 
	$(LEGION_PROF) -o ./profiler_html/$(DATE) profiling_data.gz
	mkdir -p ./profiler_html/spy_graphs/$(DATE)	
	$(LEGION_SPY) -de spy.log
	mv dataflow_main_1.pdf ./profiler_html/spy_graphs/$(DATE)
	mv event_graph_main_1.pdf ./profiler_html/spy_graphs/$(DATE)
	@echo start a webserver to view the profiling data: python -m SimpleHTTPServer
	@echo You can now find the generated pdfs under ./profiler_html/spy_graphs/
endif


# generate a profile
prof: profiling_data.gz
ifeq ($(python_version_major), 3)
	@echo Python 3 Enviroment active. Legion Profiler only supports Python 2.7 !!!
else 
	$(LEGION_PROF) -o ./profiler_html/$(DATE) profiling_data.gz
	@echo start a webserver to view the profiling data: python -m SimpleHTTPServer
endif

# generate dependency and dataflow graphs
# TODO: check python enviroment !
spy: spy.log
	mkdir -p ./profiler_html/spy_graphs/$(DATE)	
	$(LEGION_SPY) -de spy.log
	mv dataflow_main_1.pdf ./profiler_html/spy_graphs/$(DATE)
	mv event_graph_main_1.pdf ./profiler_html/spy_graphs/$(DATE)
	@echo You can now find the generated pdfs under ./profiler_html/spy_graphs/
spy-verbose: spy.log	
	mkdir -p ./profiler_html/spy_graphs/$(DATE)	
	$(LEGION_SPY) -dez spy.log
	mv dataflow_main_1.pdf ./profiler_html/spy_graphs/$(DATE)
	mv event_graph_main_1.pdf ./profiler_html/spy_graphs/$(DATE)
	@echo You can now find the generated pdfs under ./profiler_html/spy_graphs/

clean:
	rm -f spy.log
	rm -f profiling_data.gz
clean-all:
	rm -rf profiler_html
	rm -f spy.log
	rm -f profiling_data.gz