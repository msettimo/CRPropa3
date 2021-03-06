#include "crpropa/ModuleList.h"
#include "crpropa/ProgressBar.h"

#include <omp.h>
#include <algorithm>
#include <signal.h>
#ifndef sighandler_t
typedef void (*sighandler_t)(int);
#endif

using namespace std;

namespace crpropa {

bool g_cancel_signal_flag = false;
void g_cancel_signal_callback(int sig) {
	g_cancel_signal_flag = true;
}

ModuleList::ModuleList() :
		showProgress(false) {
}

ModuleList::~ModuleList() {
}

void ModuleList::setShowProgress(bool show) {
	showProgress = show;
}

void ModuleList::add(Module *module) {
	modules.push_back(module);
}

void ModuleList::process(Candidate *candidate) {
	module_list_t::iterator iEntry = modules.begin();
	while (iEntry != modules.end()) {
		ref_ptr<Module> &module = *iEntry;
		iEntry++;
		module->process(candidate);
	}
}

void ModuleList::run(Candidate *candidate, bool recursive) {
	while (candidate->isActive() && !g_cancel_signal_flag)
		process(candidate);

	// propagate secondaries
	if (recursive)
		for (size_t i = 0; i < candidate->secondaries.size(); i++)
			run(candidate->secondaries[i], recursive);
}

void ModuleList::run(candidate_vector_t &candidates, bool recursive) {
	size_t count = candidates.size();

#if _OPENMP
	std::cout << "crpropa::ModuleList: Number of Threads: " << omp_get_max_threads() << std::endl;
#endif

	ProgressBar progressbar(count);

	if (showProgress) {
		progressbar.start("Run ModuleList");
	}

	g_cancel_signal_flag = false;
	sighandler_t old_signal_handler = ::signal(SIGINT,
			g_cancel_signal_callback);

#pragma omp parallel for schedule(static, 1000)
	for (size_t i = 0; i < count; i++) {
		if (!g_cancel_signal_flag)
			run(candidates[i], recursive);

		if (showProgress)
#pragma omp critical(progressbarUpdate)
			progressbar.update();
	}

	::signal(SIGINT, old_signal_handler);
}

void ModuleList::run(Source *source, size_t count, bool recursive) {

#if _OPENMP
	std::cout << "crpropa::ModuleList: Number of Threads: " << omp_get_max_threads() << std::endl;
#endif

	ProgressBar progressbar(count);

	if (showProgress) {
		progressbar.start("Run ModuleList");
	}

	g_cancel_signal_flag = false;
	sighandler_t old_signal_handler = ::signal(SIGINT,
			g_cancel_signal_callback);

#pragma omp parallel for schedule(static, 1000)
	for (size_t i = 0; i < count; i++) {
		ref_ptr<Candidate> candidate = source->getCandidate();
		if (!g_cancel_signal_flag)
			run(candidate, recursive);

		if (showProgress)
#pragma omp critical(progressbarUpdate)
			progressbar.update();
	}

	::signal(SIGINT, old_signal_handler);
}

ModuleList::module_list_t &ModuleList::getModules() {
	return modules;
}

const ModuleList::module_list_t &ModuleList::getModules() const {
	return modules;
}

void ModuleList::showModules() const {
	crpropa::ModuleList::module_list_t::const_iterator iEntry;

	iEntry = getModules().begin();
	while (iEntry != getModules().end()) {
		const crpropa::ref_ptr<crpropa::Module> &entry = *iEntry;
		iEntry++;
		std::cout << "  " << entry->getDescription() << "\n";
	}
}

} // namespace crpropa

std::ostream &operator<<(std::ostream &out, const crpropa::ModuleList &list) {
	crpropa::ModuleList::module_list_t::const_iterator iEntry;

	iEntry = list.getModules().begin();
	while (iEntry != list.getModules().end()) {
		const crpropa::ref_ptr<crpropa::Module> &entry = *iEntry;
		iEntry++;
		out << "  " << entry->getDescription() << "\n";
	}
	return out;
}
