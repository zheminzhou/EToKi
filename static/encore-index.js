var formatCount = d3.format(".0f");

function InteractPlot(div) {
	var self = this;
	this.container = div;
	div.addClass('framed')
		.append('<div id="menu" class="tab-handle" style="display:none;position:absolute"></div>')
		.on('click', function(d) {
			self.menu.hide(300);
		});
	this.menu = self.container.find("#menu")
	.css('width', div.width())
	.append('<button id="button-default" class="tab-item-handle" value="default">Default</button>');
	this.menu.find("#button-default").on('click', function(e, ui) {
		self.plotDim.x = self.plotDim.x0;
		self.plotDim.y = self.plotDim.y0;
		self.plotDim.dw = self.plotDim.dh = 1;
		self.draw();
	})
	this.container.on("contextmenu", function(e, ui) {
		e.preventDefault();
		self.menu.show(300);
	})
	this.inSelection = false;
	this.margin = 50;
	this.plotDim = {
		x0 : this.margin, 
		y0 : div.height() - this.margin, 
		x1 : div.width() - this.margin,
		y1 : this.margin,
		x : this.margin, 
		y : div.height() - this.margin, 
		w : div.width() - 2*this.margin, 
		h : 2*this.margin-div.height(), 
		dw : 1,
		dh : 1,
	}
	div.append('<svg id="plotView"></svg>');
	this.plotView = d3.select(div.find('#plotView').selector)
		.attr("height", div.height())
		.attr("width", div.width())
		.call(d3.drag().on('start', function(d) {
				if (d3.event.sourceEvent.shiftKey) {
					self.drawSelection(d3.event.sourceEvent.offsetX, d3.event.sourceEvent.offsetY, null, null);
				}
			})
			.on('drag', function(d) {
				if (! self.inSelection) {
					self.plotDim.x = self.plotDim.x+d3.event.dx
					self.plotDim.y = self.plotDim.y+d3.event.dy
					self.draw();
				} else {
					self.drawSelection(null, null, d3.event.sourceEvent.offsetX, d3.event.sourceEvent.offsetY);
				}
			})
			.on('end', function(d) {
				self.drawSelection(null, null, null, null);
			})
		)
		.call(d3.zoom().on('zoom', function(d) {
			if (! d3.event.sourceEvent) {
				self.mat.forEach(function(d) {d.style.selected=false;});
				self.plotView.selectAll('.selected').remove();
				return;
			}
			var scale = d3.event.sourceEvent.deltaY < 0 ? 1.1 : 1/1.1;
			if (d3.event.sourceEvent.shiftKey) {
				self.plotDim.dh *= scale;
				self.plotDim.y = d3.event.sourceEvent.offsetY - (d3.event.sourceEvent.offsetY-self.plotDim.y)*scale;
			} else if (d3.event.sourceEvent.ctrlKey) {
				self.plotDim.dh *= scale;
				self.plotDim.y = d3.event.sourceEvent.offsetY - (d3.event.sourceEvent.offsetY-self.plotDim.y)*scale;
				self.plotDim.dw *= scale;
				self.plotDim.x = d3.event.sourceEvent.offsetX - (d3.event.sourceEvent.offsetX-self.plotDim.x)*scale;
			} else {
				self.plotDim.dw *= scale;
				self.plotDim.x = d3.event.sourceEvent.offsetX - (d3.event.sourceEvent.offsetX-self.plotDim.x)*scale;
			}
			self.draw();
		}));
};

InteractPlot.prototype.drawSelection = function(x1, y1, x2, y2) {
	var self = this;
	if ( x2 === null ) {
		if (x1 === null) {
			self.inSelection = false;
			var selectRegion = self.plotView.selectAll('.selection');
			if (selectRegion.size() == 0) {return;}
			var selection = selectRegion.data()[0];
			selectRegion.remove();
			if (selection[0] > selection[2]) {[selection[0], selection[2]] = [selection[2], selection[0]]};
			if (selection[1] > selection[3]) {[selection[1], selection[3]] = [selection[3], selection[1]]};
			var selected = self.plotView.selectAll('.data').filter(function(d) {
				if (d.r === undefined) {
					var dim = [d.x, d.y, d.x+d.w, d.y+d.h];
					var ovlX = Math.min(dim[2], selection[2]) - Math.max(dim[0], selection[0]);
					if (ovlX >= 0) {
						var ovlY = Math.min(dim[3], selection[3]) - Math.max(dim[1], selection[1]);
						if (ovlY >= 0) {
							return true;
						}
					}
				} else {
					if (d.x >= selection[0] && d.x <= selection[2]) {
						if (d.y >= selection[1]-d.r && d.y <= selection[3]+d.r) {
							return true;
						}
					} else if (d.y >= selection[1] && d.y <= selection[3]) {
						if (d.x >= selection[0]-d.r && d.x <= selection[2]+d.r) {
							return true;
						}
					} else {
						var mx = Math.min(Math.abs(selection[0]-d.x), Math.abs(selection[2]-d.x));
						var my = Math.min(Math.abs(selection[1]-d.y), Math.abs(selection[3]-d.y));
						if (mx <= d.r && my <=d.r && d.r*d.r >= mx*mx + my*my ) {
							return true;
						}
					}
				}
				return false;
			});
			if (selected.size() > 0) {
				if (selected.filter(function(d) {return d.style.selected != true;}).size() > 0) {
					selected.each(function(d) {d.style.selected = true;});
				} else {
					selected.each(function(d) {d.style.selected = false;});
				}
				self.draw();
				var datum = []
				self.mat.forEach(function(d) {
					if (d.style.selected === true) {
						datum.push(d.data);
					};
				});
				$(document).trigger("plot:selected", {datum:datum});
			}
		} else {
			self.inSelection = true;
			var selection = [[x1, y1, x1, y1]];
			self.plotView.selectAll('.selection')
				.data(selection)
				.enter().append('rect')
					.attr('class', 'selection')
					.attr('x', function(d){return d[0];})
					.attr('y', function(d){return d[1];})
					.attr('width', function(d){return d[2]-d[0];})
					.attr('height', function(d){return d[3]-d[1];})
					.attr('fill', 'green')
					.attr('opacity', '0.5')
		}
	} else {
		var selectRegion = self.plotView.selectAll('.selection');
		var selection = selectRegion.data()[0];
		selection[2] = x2;
		selection[3] = y2;
		if (selection[2] > selection[0]) {
			selectRegion.attr('width', function(d){return d[2]-d[0];});
		} else {
			selectRegion.attr('x', function(d){return d[2];})
				.attr('width', function(d){return d[0]-d[2];});
		}
		if (selection[3] > selection[1]) {
			selectRegion.attr('height', function(d){return d[3]-d[1];});
		} else {
			selectRegion.attr('y', function(d){return d[3];})
				.attr('height', function(d){return d[1]-d[3];});
		}
	}
}

InteractPlot.prototype.draw = function(option) {
	option = option ? option : {};
	var data   = option.data     ?   option.data : null;
	var types  = option.type     ?   option.type.split(',') : (this.types ? this.types : ['histogram']);
	var x      = option.x        ?      option.x : function(d) {return d[0]?d[0]:1;};
	var y      = option.y        ?      option.y : function(d) {return d[1]?d[1]:1;};
	var z      = option.z        ?      option.z : function(d) {return 1;};
	var style  = option.style    ?  option.style : function(d) {return {stroke:'black', 'stroke-width':0.5, 'fill':'steelblue'}; }

	var self = this;
	this.plotView.selectAll('*').remove();

	if (data !== null) {
		var mat = data.map(function(d) {
			return {x:x(d), y:y(d), z:z(d), style:style(d), data:d};
		});
		self.mat = mat;
	} else {
		var mat = self.mat;
	}
	var style_list = {};
	mat.forEach(function(d) { 
		Object.keys(d.style).forEach(function(k) {
			style_list[k] = 1;
		}); 
	});
	self.types = types;
	var assign_style = function(item, style_list) {
		Object.keys(style_list).forEach(function(k) {
			item.attr(k, function(d) {return d.style[k];});
		})
	}

	var xx = d3.scaleLinear().domain([0, 1.05*Math.max.apply(
		Math, mat.map(function(d) {return d.x})
	)]).range([self.plotDim.x, self.plotDim.x+self.plotDim.w*self.plotDim.dw]);
	var yy = d3.scaleLinear().domain([0, 1.05*Math.max.apply(
		Math, mat.map(function(d) {return d.y})
	)]).range([self.plotDim.y, self.plotDim.y+self.plotDim.h*self.plotDim.dh]);

	for (var id in types) {
		type = types[id];
		if (type === 'histogram') {
			var bar = this.plotView.selectAll('.data')
				.data(mat.map(function(d) {return {
					x : xx(d.x),
					y : yy(d.y),
					w : xx(d.z) - xx(0), 
					h : self.plotDim.y - yy(d.y),
					style : d.style, 
					data  : d.data, 
				};}))
				.enter().append('g')
					.attr("class", "data");
			var bar2 = bar.append('rect')
					.attr("x", function(d) {return d.x;})
					.attr("y", function(d) {return d.y;})
					.attr("width", function(d) {return d.w;})
					.attr("height", function(d) {return d.h;});
			assign_style(bar2, style_list);
			bar.filter(function(d) {return d.style.selected === true})
				.append('rect')
				.attr("class", "selected")
				.attr("x", function(d) {return d.x;})
				.attr("y", function(d) {return d.y;})
				.attr("width", function(d) {return d.w;})
				.attr("height", function(d) {return d.h;})
				.attr("opacity", 0.4)
				.attr("fill", "red")
				.attr("stroke-width", 3)
				.attr("stroke", "black");
		} else if (type === 'scatter') {
			mat.forEach(function(d) {
				d[type] = {
					x : xx(d.x), y: yy(d.y), r: 3*d.z, 
					w : 0, h: 0,
				};
			});

			var dot = this.plotView.selectAll('.data')
				.data(mat.map(function(d) {return {
					x : xx(d.x),
					y : yy(d.y),
					r: 3*d.z, 
					style : d.style, 
					data  : d.data, 
				};}))
				.enter().append('g')
					.attr('class', 'data');
					
			var dot2 = dot.append('circle')
					.attr('cx', function(d) {return d.x;})
					.attr('cy', function(d) {return d.y;})
					.attr('r', function(d) {return d.r;});
			assign_style(dot2, style_list);
			dot.filter(function(d) {return d.style.selected === true})
				.append('circle')
				.attr("class", "selected")
				.attr('cx', function(d) {return d.x;})
				.attr('cy', function(d) {return d.y;})
				.attr('r', function(d) {return d.r;})
				.attr("opacity", 0.4)
				.attr("fill", "red")
				.attr("stroke-width", 3)
				.attr("stroke", "black");
		} else if (type === 'line') {
			this.plotView.append("path")
			  .datum(mat)
			  .attr(d.style)
			  .attr("d", d3.line()
							.x(function(d) {return xx(d.x);})
							.y(function(d) {return yy(d.y);})
				);
		}
	}
	this.plotView.append('g')
		.attr('transform', "translate(0, "+self.plotDim.y0+")")
		.call(d3.axisBottom(xx))
	this.plotView.append('g')
		.attr('transform', "translate("+self.plotDim.x0+",0)")
		.call(d3.axisLeft(yy))
	this.plotView.append('g')
		.attr('transform', "translate(0, "+self.plotDim.y1+")")
		.call(d3.axisTop(xx))
	this.plotView.append('g')
		.attr('transform', "translate("+self.plotDim.x1+",0)")
		.call(d3.axisRight(yy))

	return this;
};

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
	}) );
	var presencePlot = new InteractPlot( $("#presence").css({
		height: 500, width: 900
	}) );
	var pseudoPlot = new InteractPlot( $("#cds").css({
		height: 500, width: 900
	}) );
	var variPlot = new InteractPlot( $("#variation").css({
		height: 500, width: 900
	}) );

	$(document).on("redraw", function(e, data) {
		genomePlot.draw( {
			data : function(d){
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
			x : function(d) {return d[0]}, 
			y : function(d) {return d[1].length},
			type : 'histogram'
		});
		presencePlot.draw( {
			data : function(d){
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
			x : function(d) {return d[0]}, 
			y : function(d) {return d[1].length},
			type : 'histogram'
		});
		pseudoPlot.draw( {
			data : function(d){
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
			x : function(d) {return d[0]}, 
			y : function(d) {return d[1].length},
			type : 'histogram'
		});

		variPlot.draw( {
			data : data.locusStatus.filter(function(v) {return v.Activated;}), 
			x : function(d) {return Number(d.n_Size)}, 
			y : function(d) {return Number(d.n_Allele)},
			type : 'scatter'
		});
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
