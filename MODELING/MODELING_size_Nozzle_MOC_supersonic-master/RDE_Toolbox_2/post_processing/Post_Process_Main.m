function Post_Process_Main(PostProc)
%POST_PROCESS_MAIN Entry point for all post-processing / plotting steps.
%   Post_Process_Main(PostProc) takes a struct PostProc holding whatever
%   each step needs, and runs the post-processing steps in order. Called
%   from MoC_Main.m once the flow-field computation is finished.
%
%   PostProc fields used so far:
%       xMax - domain extent in x
%       yMax - domain extent in y

ax = Setup_Graph_Domain(PostProc.xMax, PostProc.yMax);

% Future stuff goes here
end