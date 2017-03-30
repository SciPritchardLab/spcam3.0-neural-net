!  CRM domain grid-size:

	  integer YES3D
	  integer crm_nx, crm_ny, crm_nz
	  real crm_dx, crm_dy, crm_dt
          parameter (YES3D = 0)   ! 1 - 3D CRM, 0 - 2D CRM
	  parameter (crm_nx = 32, crm_ny = 1, crm_nz = 28)
	  parameter (crm_dx = 4000., crm_dy = crm_dx)
	  parameter (crm_dt = 20.)
