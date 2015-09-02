function params_final = generate_ba_instance_params
% generate_ba_instance_params:
%       assuming 11 camera parameters and zach robustifier
%       ntasks x [n_cams, n_pts, n_obs, total_in_params]
%
%       returns a cell array of:
%           [num_cams, num_pts, num_observations, num_in_parameters,...
%               num_out_parameters, ]

ncam_params = 11;

params = {
          [21    11315   36455    ],...
          [49    7776    31843    ],...
          [88    64298   383937   ],...
          [93    61203   287451   ],...
          [161   48126   182072   ],...
          [245   198739  1091386  ],...
          [253   163691  899155   ],...
          [257   65132   225911   ],...
          [356   226730  1255268  ],...
          [372   47423   204472   ],...
          [810   88814   393775   ],...
          [539   65220   277273   ],...
          [1102  780462  4052340  ],...
          [1197  126327  563734   ],...
          [1544  942409  4750193  ],...
          [1723  156502  678718   ],...
          [1778  993923  5001946  ],...
          [1936  649673  5213733  ],...
          [4585  1324582 9125125  ],...
          [13682 4456117 28987644 ]...
          };

p = reshape(cell2mat(params),[3 numel(params)])';

% cams, points and weights
nin = ncam_params*p(:,1) + 3*p(:,2) + p(:,3);
% reprojection error, weights
nout = 2*p(:,3) + p(:,3);

% nonzero elements (principal point contributes by one coordinate to one
% row)
nnonzero = (ncam_params-1 + 3 + 1) * p(:,3) + p(:,3);

[~,order]=sort(nnonzero);
params_final = {};
for i=order'
    params_final{i} = [params{i} nin(i) nout(i) nnonzero(i)];
end

end