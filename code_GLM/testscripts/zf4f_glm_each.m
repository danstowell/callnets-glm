% script to run the GLM fitting on a sequence of data chunks from CSV files of calling data

plotcols = {'r', 'b', 'g', 'm'};

% CONFIG: paths to source-data files and output data files. (the output data files are also used as source when we re-analyse a resimulated data series.)
zf4fdir = '../data/zf4f'
gilldir = '../data/gilletal_trial2'
csvoutdir = 'outcsv'

setses = { ...
struct('runname', 'session2fulla', 'indexmapper', [4,3,2,1], 'startsecs',  300, 'endsecs', 1200, 'resimuldur', 0, 'datapath', sprintf('%s/zcompiled_session2full', zf4fdir)), ...
struct('runname', 'session2fullb', 'indexmapper', [4,3,2,1], 'startsecs', 1200, 'endsecs', 2100, 'resimuldur', 0, 'datapath', sprintf('%s/zcompiled_session2full', zf4fdir)), ...
struct('runname', 'session2fullc', 'indexmapper', [4,3,2,1], 'startsecs', 2100, 'endsecs', 3000, 'resimuldur', 0, 'datapath', sprintf('%s/zcompiled_session2full', zf4fdir)), ...
struct('runname', 'session2fulld', 'indexmapper', [4,3,2,1], 'startsecs', 3000, 'endsecs', 3900, 'resimuldur', 0, 'datapath', sprintf('%s/zcompiled_session2full', zf4fdir)), ...
struct('runname', 'session3fulla', 'indexmapper', [1,2,3,4], 'startsecs',  350, 'endsecs', 1250, 'resimuldur', 0, 'datapath', sprintf('%s/zcompiled_session3full', zf4fdir)), ...
struct('runname', 'session3fullb', 'indexmapper', [1,2,3,4], 'startsecs', 1250, 'endsecs', 2150, 'resimuldur', 0, 'datapath', sprintf('%s/zcompiled_session3full', zf4fdir)), ...
struct('runname', 'session3fullc', 'indexmapper', [1,2,3,4], 'startsecs', 2150, 'endsecs', 3050, 'resimuldur', 0, 'datapath', sprintf('%s/zcompiled_session3full', zf4fdir)), ...
struct('runname', 'session3fulld', 'indexmapper', [1,2,3,4], 'startsecs', 3050, 'endsecs', 3950, 'resimuldur', 0, 'datapath', sprintf('%s/zcompiled_session3full', zf4fdir)), ...
struct('runname', 'session2full', 'indexmapper', [4,3,2,1], 'startsecs',  300, 'endsecs', 3900, 'resimuldur', 600, 'datapath', sprintf('%s/zcompiled_session2full', zf4fdir)), ...
struct('runname', 'session3full', 'indexmapper', [1,2,3,4], 'startsecs',  350, 'endsecs', 3950, 'resimuldur', 600, 'datapath', sprintf('%s/zcompiled_session3full', zf4fdir)), ...
%
struct('runname', 'session2full_resim', 'indexmapper', [4,3,2,1], 'startsecs',  0, 'endsecs', 60000, 'resimuldur', 0, 'datapath', sprintf('%s/zf4f_glm_stats_session2fullsof_resimulated', csvoutdir)), ...
struct('runname', 'session3full_resim', 'indexmapper', [4,3,2,1], 'startsecs',  0, 'endsecs', 60000, 'resimuldur', 0, 'datapath', sprintf('%s/zf4f_glm_stats_session3fullsof_resimulated', csvoutdir)), ...
%
struct('runname', 'gill1ind', 'indexmapper', 1:8, 'startsecs', 0, 'endsecs', 14400, 'resimuldur', 0, 'datapath', sprintf('%s/gill_rawfile_1ind', gilldir)), ...
struct('runname', 'gill4ind', 'indexmapper', 1:8, 'startsecs', 0, 'endsecs', 14400, 'resimuldur', 0, 'datapath', sprintf('%s/gill_rawfile_4ind', gilldir)), ...
struct('runname', 'gill5ind', 'indexmapper', 1:8, 'startsecs', 0, 'endsecs', 14400, 'resimuldur', 0, 'datapath', sprintf('%s/gill_rawfile_5ind', gilldir)), ...
struct('runname', 'gill6ind', 'indexmapper', 1:8, 'startsecs', 0, 'endsecs', 14400, 'resimuldur', 0, 'datapath', sprintf('%s/gill_rawfile_6ind', gilldir)), ...
%
struct('runname', 'gill1type', 'indexmapper', 1:40, 'startsecs', 0, 'endsecs', 14400, 'resimuldur', 0, 'datapath', sprintf('%s/gill_rawfile_1type', gilldir)), ...
struct('runname', 'gill4type', 'indexmapper', 1:40, 'startsecs', 0, 'endsecs', 14400, 'resimuldur', 0, 'datapath', sprintf('%s/gill_rawfile_4type', gilldir)), ...
struct('runname', 'gill5type', 'indexmapper', 1:40, 'startsecs', 0, 'endsecs', 14400, 'resimuldur', 0, 'datapath', sprintf('%s/gill_rawfile_5type', gilldir)), ...
struct('runname', 'gill6type', 'indexmapper', 1:40, 'startsecs', 0, 'endsecs', 14400, 'resimuldur', 0, 'datapath', sprintf('%s/gill_rawfile_6type', gilldir)), ...
};


nlfuns = {@softplus, @expfun};
for whichnlf = 1:length(nlfuns)
	nlfun = nlfuns{whichnlf};
	nlname = func2str(nlfun)(1:3);
	numcalls   = struct();
	resultspos = struct();
	resultsval = struct();
	negloglis  = struct();
	dcs        = struct();
	for whichset=1:size(setses,2)
		d = setses{whichset};
		if whichnlf==1 || whichnlf==2
			resimuldur = d.resimuldur;
		else
			resimuldur = 0;
		end
		runlabel = sprintf('%s%s', d.runname, func2str(nlfun)(1:3));
		csvpath = sprintf('%s.csv', d.datapath);

		disp(sprintf('Fitting with nonlin %s on %s', func2str(nlfun), csvpath));

		k = length(d.indexmapper);
		regln = -1; % NOTE default regularisation strength here
		[numcalls.(d.runname), resultspos.(d.runname), resultsval.(d.runname), negloglis.(d.runname), dcs.(d.runname)] = dofit_fromcsv_GLM_zf4f(csvpath, runlabel, k, d.indexmapper, d.startsecs, d.endsecs, regln, 'outplot', csvoutdir, d.resimuldur, nlfun);
	end
end;

