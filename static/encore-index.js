var formatCount = d3.format(".0f");

function homeTab(div, option) {
	var self = this;
	option = option ? option : {};
	option.autoHidden = option.autoHidden ? option.autoHidden : true;
	option.position   =   option.position ? option.position   : 'right';
	option.size       =       option.size ? option.size       : 100;
	
	this.div = div;
	if (option.position)
}

function floatingTab(div, tabs) {
	var self = this;
	this.div = div
		.addClass('floating-tab')
		.append('<div class="tab-handle"></div>')

	this.handle = this.div.find(".tab-handle");

	this.dim = {width:0, height:0};
	for (var i in tabs) {
		[name, tab] = tabs[i];
		this.div.append(tab.addClass('tab-item'));
		if (tab.width() > this.dim.width) {
			this.dim.width = tab.width();
		}
		if (tab.height() > this.dim.height) {
			this.dim.height = tab.height();
		}
		tab.hide();
		this.handle.append('<button class="tab-item-handle" value="'+tab.selector+'">'+name+'</button>');
	}
	this.div.css('width', this.dim.width);
	this.handle.on('click', function(e, ui) {
		self.div.find(".tab-item").hide(300);
		self.div.find(".tab-item-handle").removeClass('active');
		self.div.animate({height: 30}, 300);
	});
	this.handle.find(".tab-item-handle")
		.on("click", function(e, ui) {
			e.stopPropagation();
			var tab = $($(this).attr("value"));
			if (tab.css('display') === 'none') {
				self.div.find(".tab-item-handle").removeClass('active');
				$(this).addClass('active');
				self.div.animate({'height': self.dim.height+30}, 300).find(".tab-item").hide();
				$($(this).attr("value")).show();
			}
		});
}

function nGrid(div, cols) {
	var self = this;
	this.filtedData = [];
  	this.columns = cols.map(function(c) {
		return {id:c, name:c, field:c, width: 100, sortable:true, editor: Slick.Editors.Text};
	});;
  	this.options = {
  		autoEdit: false,
  		enableCellNavigation: true,
  		topPanelHeight: 50,
	};
  	this.grid = new Slick.Grid(div, [], this.columns, this.options);
  	this.grid.setSelectionModel(new Slick.RowSelectionModel());
	div.append('<ul id="contextMenu" style="display:none;position:absolute"><b>:</b><li data="activate">Activate</li><li data="inactivate">Inactivate</li></ul>');
	div.find("#contextMenu").click(function (e) {
		if (!$(e.target).is("li")) {
		  return;
		}
		var rows = $(this).data("rows").map(function(i) {return self.grid.getData()[i].id});
		$(document).trigger($(e.target).attr('data'), {ids:rows});
	});

	this.grid.onContextMenu.subscribe(function (e, ui) {
		e.preventDefault();
		var cur_row = ui.grid.getCellFromEvent(e).row;
		var selected_rows = ui.grid.getSelectedRows();
		if (selected_rows.filter(function(r) {return r === cur_row;}).length === 0) {
			ui.grid.setSelectedRows([cur_row]);
			selected_rows = [cur_row];
		}

		div.find("#contextMenu")
			.data('rows', selected_rows)
			.css("top", e.pageY - div.offset().top)
			.css("left", e.pageX - div.offset().left)
			.show();
		$(document).on("click", function () {
			div.find("#contextMenu").hide();
		});
	});

	this.grid.onSort.subscribe(function (e, args) {
		var currentSortCol = args.sortCol.id;
		var isAsc = args.sortAsc ? 1 : -1;

		self.grid.getData().sort(function(d1, d2) {
			var [v1, v2] = [d1[currentSortCol], d2[currentSortCol]];
			var [t1, t2] = [isNaN(v1), isNaN(v2)];
			if (t1 === t2) {
				if (t1 === true) {
					return v1 === v2 ? (d1.__index - d2.__index) : (v1 > v2 ? isAsc : -isAsc) ;
				} else {
					var [v1, v2] = [Number(v1), Number(v2)];
					return v1 === v2 ? (d1.__index - d2.__index) : isAsc*(v1 - v2) ;
				}
			} else {
				return isAsc*(t1 > t2);
			}
		});

		var t = 0;
		self.grid.getData().forEach(function (d) {t = self.grid.idOrder[d.id] = t+1;});
		self.grid.invalidateAllRows();
		self.grid.render();
	});
}

nGrid.prototype.redraw = function(data) {
	var self = this;
	if (! self.grid.idOrder) {
		var t=0;
		this.grid.idOrder = {};
		data.forEach(function (d) {t = self.grid.idOrder[d.id] = t+1;});
	} else {
		data.sort(function(d1, d2) {
			return self.grid.idOrder[d1.id] - self.grid.idOrder[d2.id];
		});
	}
	this.grid.setData(data);
	this.grid.invalidateAllRows();
	this.grid.render();
};

function centralMatrix(data) {
	var self = this;
	this.g = JSON.parse(JSON.stringify(data));
	this.l = {};
	var gg = this.g,
		ll = this.l;

	this.active = {};
	Object.keys(self.g).map(function(g) {
		self.active[g] = 1;
		Object.keys(self.g[g]).map(function(l) {
			self.active[l] = 1;
			self.g[g][l][3] = self.g[g][l][0];
			self.l[l] = {};
		});
	});
	Object.keys(self.g).map(function(g) {
		Object.keys(self.g[g]).map(function(l) {
			self.l[l][g] = self.g[g][l];
		})
	});
	this.GenomeStatus();
	this.LocusStatus();
	$(document).trigger("redraw", this);

	$(document).on("activate", function(e, ui) {
		self.unfilt(ui.ids);
		self.GenomeStatus();
		self.LocusStatus();
		$(document).trigger("redraw", self);
	})
	.on("inactivate", function(e, ui) {
		self.filt(ui.ids);
		self.GenomeStatus();
		self.LocusStatus();
		$(document).trigger("redraw", self);
	})

};

centralMatrix.prototype.filt = function(idList) {
	var self = this;
	idList.map(function(i) {
		self.active[i] = 0;
	});
}

centralMatrix.prototype.unfilt = function(idList) {
	var self = this;
	idList.map(function(i) {
		self.active[i] = 1;
	});
}
centralMatrix.prototype.GenomeStatus = function(idList) {
	if (idList === undefined) {
		idList = Object.keys(this.g);
	}
	var self = this;
	this.genomeStatus = idList.map(function(g) {
		var v = Object.keys(self.g[g]).filter(function(l) {return self.active[l] > 0;}).map(function(l) {return self.g[g][l];});
		nTotal = v.length;
		var v = v.filter(function(d) {return d[0]>0;});
		nPresence = v.length;
		nBase = v.map(function(d) {return d[1];}).reduce(function(v1, v2) {return v1+v2;});
		nIntact = v.map(function(d) {return d[0]>2;}).reduce(function(v1, v2) {return v1+v2;});
		return {Activated:self.active[g], id:g, n_Total:nTotal, n_Presence:formatCount(100*nPresence/nTotal), n_Intact:formatCount(100*nIntact/nPresence), n_Base:nBase};
	});
	return this.genomeStatus;
}
centralMatrix.prototype.LocusStatus = function(idList) {
	if (idList === undefined) {
		idList = Object.keys(this.l);
	}
	var self = this;
	this.locusStatus = idList.map(function(l) {
		var v = Object.keys(self.l[l]).filter(function(g) {return self.active[g] > 0;}).map(function(g) {return self.l[l][g];});
		nTotal = v.length;
		var v = v.filter(function(d) {return d[0]>0;});
		nPresence = v.length;
		allele_ids = {};
		v.forEach(function(d) {allele_ids[d[2]] = 1;});
		nAllele = Object.keys(allele_ids).length;
		nBase = v.length ? v.map(function(d) {return d[1];}).reduce(function(v1, v2) {return v1+v2;}) : 0;
		nIntact = v.length ? v.map(function(d) {return d[0]>2;}).reduce(function(v1, v2) {return v1+v2;}) : 0;
		return {Activated:self.active[l], id:l, n_Total:nTotal, n_Presence:formatCount(100*nPresence/nTotal), n_Intact:formatCount(100*nIntact/nPresence), n_Allele:nAllele, n_Size:formatCount(nBase/nPresence)};
	});
	return this.locusStatus;
}


$( function() {
	var genomePlot = new InteractPlot( $("#genome").css({
		height: 500, width: 900
	}).resizable().draggable() );
	var presencePlot = new InteractPlot( $("#presence").css({
		height: 500, width: 900
	}).resizable().draggable() );
	var pseudoPlot = new InteractPlot( $("#cds").css({
		height: 500, width: 900
	}).resizable().draggable() );
	var variPlot = new InteractPlot( $("#variation").css({
		height: 500, width: 900
	}).resizable().draggable() );

	$(document).on("redraw", function(e, data) {
		genomePlot.histogram( 
			function(d){
				var p = {};
				d.filter(function(v) {return v.Activated;}).forEach(function(v) {
					var prop = v.n_Presence;
					if ( p[prop] ) {
						p[prop].push(v);
					} else {
						p[prop] = [v];
					}
				});
				return Object.keys(p).map(function(prop) {return [Number(prop), p[prop]]});
			}(data.genomeStatus), 
			{
				x : function(d) {return d[0]}, 
				y : function(d) {return d[1].length},
		}).draw();
		presencePlot.histogram(
			function(d){
				var p = {};
				d.filter(function(v) {return v.Activated;}).forEach(function(v) {
					var prop = v.n_Presence;
					if ( p[prop] ) {
						p[prop].push(v);
					} else {
						p[prop] = [v];
					}
				});
				return Object.keys(p).map(function(prop) {return [Number(prop), p[prop]]});
			}(data.locusStatus), 
			{
				x : function(d) {return d[0]}, 
				y : function(d) {return d[1].length},
				type : 'histogram'
		}).draw();
		pseudoPlot.histogram(
			function(d){
				var p = {};
				d.filter(function(v) {return v.Activated;}).forEach(function(v) {
					var prop = v.n_Intact;
					if ( p[prop] ) {
						p[prop].push(v);
					} else {
						p[prop] = [v];
					}
				});
				return Object.keys(p).map(function(prop) {return [Number(prop), p[prop]]});
			}(data.locusStatus), 
			{
				x : function(d) {return d[0]}, 
				y : function(d) {return d[1].length},
		}).draw();

		variPlot.scatter(
			data.locusStatus.filter(function(v) {return v.Activated;}), 
			{
				x : function(d) {return Number(d.n_Size)}, 
				y : function(d) {return Number(d.n_Allele)},
				r : function(d) {return 10},
		}).draw();
		if (locusGrid) {
			locusGrid.redraw(data.locusStatus);
		};
		if (genomeGrid) {
			genomeGrid.redraw(data.genomeStatus);
		}

	});

	var genomeGrid = new nGrid( $("#genome-grid").css({
		height: 600, width: 650,
	}), ['Activated', 'id', 'n_Total', 'n_Base', 'n_Presence', 'n_Intact']);
	var locusGrid = new nGrid( $("#locus-grid").css({
		height: 600, width: 650,
	}), ['Activated', 'id', 'n_Total', 'n_Allele', 'n_Size', 'n_Presence', 'n_Intact']);

	data = new centralMatrix(matrix);

	gridTab = new floatingTab( $("#grid"), [
		["Genomes", $("#genome-grid")], 
		["Loci", $("#locus-grid")], 
	] );

});
