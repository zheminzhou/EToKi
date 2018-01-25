function InteractPlot(div) {
	var self = this;
	this.container = div;
	
	this.inSelection = false;
	this.margin = 50;
	this.plotDim = {
		x0 : this.margin, 
		y0 : this.container.height() - this.margin, 
		x1 : this.container.width() - this.margin,
		y1 : this.margin,
		x  : this.margin, 
		y  : this.container.height() - this.margin, 
		w  : this.container.width() - 2*this.margin, 
		h  : 2*this.margin - this.container.height(), 
		dw : 1,
		dh : 1,
	}
	this.addMenu();
	this.container.append('<svg id="plotView"></svg>');
	this.plotView = d3.select(this.container.find('#plotView').selector)
		.attr("height", this.container.height())
		.attr("width", this.container.width())
		.call(d3.drag().on('start', function(d) {
				if (d3.event.sourceEvent.shiftKey) {
					self.dragSelect('start', d3.event.sourceEvent.offsetX, d3.event.sourceEvent.offsetY);
				}
			})
			.on('drag', function(d) {
				if (! self.inSelection) {
					self.plotDim.x = self.plotDim.x+d3.event.dx
					self.plotDim.y = self.plotDim.y+d3.event.dy
					self.draw();
				} else {
					self.dragSelect('drag', d3.event.sourceEvent.offsetX, d3.event.sourceEvent.offsetY);
				}
			})
			.on('end', function(d) {
				self.dragSelect('end');
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

InteractPlot.prototype.addMenu = function() {
	var self = this;
	this.container.addClass('framed')
		.append('<div id="menu" class="tab-handle" style="display:none;position:absolute"></div>')
		.on("contextmenu", function(e, ui) {
			e.preventDefault();
			self.menu.show(300);
		})
		.on('click', function(d) {
			self.menu.hide(300);
		});
	this.menu = this.container.find("#menu")
		.append('<button id="button-default" class="tab-item-handle" value="default">Default</button>');
	this.menu.find("#button-default").on('click', function(e, ui) {
		self.plotDim.x = self.plotDim.x0;
		self.plotDim.y = self.plotDim.y0;
		self.plotDim.dw = self.plotDim.dh = 1;
		self.draw();
	});
};


InteractPlot.prototype.dragSelect = function(event, x, y) {
	var self = this;
	if (event == 'end') {
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
	} else if (event == 'start') {
		self.inSelection = true;
		var selection = [[x, y, x, y]];
		self.plotView.selectAll('.selection')
			.data(selection)
			.enter().append('rect')
				.attr('class', 'selection')
				.attr('x', function(d){return d[0];})
				.attr('y', function(d){return d[1];})
				.attr('width', function(d){return d[2]-d[0];})
				.attr('height', function(d){return d[3]-d[1];})
				.attr('fill', 'green')
				.attr('opacity', '0.5');
	} else if (event == 'drag') {
		var selectRegion = self.plotView.selectAll('.selection');
		var selection = selectRegion.data()[0];
		selection[2] = x;
		selection[3] = y;
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
	var x      = option.x        ?      option.x : function(d) {return d[0] ? d[0] : 1;};
	var y      = option.y        ?      option.y : function(d) {return d[1] ? d[1] : 1;};
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
	self.types = types;
	
	var assign_style = function(item, style_list) {
		var style_list = {};
		item.each(function(d) {
			Object.keys(d.style).forEach( {
				function(k) {
					style_list[k] = 1;
				}
			});
		});
		Object.keys(style_list).forEach(function(k) {
			item.attr(k, function(d) {return d.style[k];});
		})
	}

	var plot_selected = function(item) {
		if (item.data()[0].r) {
			var drawing = item.append('circle')
					.attr('cx', function(d) {return d.x;})
					.attr('cy', function(d) {return d.y;})
					.attr('r', function(d) {return d.r;});

			item.filter(function(d) {return d.style.selected === true})
				.append('circle')
				.attr("class", "selected")
				.attr('cx', function(d) {return d.x;})
				.attr('cy', function(d) {return d.y;})
				.attr('r', function(d) {return d.r;})
				.attr("opacity", 0.4)
				.attr("fill", "red")
				.attr("stroke-width", 3)
				.attr("stroke", "black");
		} else {
			var drawing = item.append('rect')
					.attr("x", function(d) {return d.x;})
					.attr("y", function(d) {return d.y;})
					.attr("width", function(d) {return d.w;})
					.attr("height", function(d) {return d.h;});

			item.filter(function(d) {return d.style.selected === true})
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
		}
		assign_style(drawing);
	}
	var xx = d3.scaleLinear().domain([0, 1.05*Math.max.apply(
		Math, mat.map(function(d) {return d.x})
	)]).range([self.plotDim.x, self.plotDim.x+self.plotDim.w*self.plotDim.dw]);
	var yy = d3.scaleLinear().domain([0, 1.05*Math.max.apply(
		Math, mat.map(function(d) {return d.y})
	)]).range([self.plotDim.y, self.plotDim.y+self.plotDim.h*self.plotDim.dh]);

	for (var id in types) {
		type = types[id];
		if (type !== 'line') {
			if (type === 'histogram') {
				var item = this.plotView.selectAll('.data')
					.data(mat.map(function(d) {return {
							x : xx(d.x),
							y : yy(d.y),
							w : xx(d.z) - xx(0), 
							h : self.plotDim.y - yy(d.y),
							style : d.style, 
							data  : d.data, 
						};
					}))
					.enter().append('g')
						.attr("class", "data");
			} else if (type === 'stack') {
				var y_stack = {};
				for (var i in mat) {
					var d = mat[i];
					if (! y_stack[d.x]) {
						d.y1 = 0;
						y_stack[d.x] = d.y;
					} else {
						d.y1 = y_stack[d.x];
						[y_stack[d.x], d.y] = [d.y, d.y + d.y1];
					}
				};
				var item = this.plotView.selectAll('.data')
					.data(mat.map(function(d) {return {
							x : xx(d.x),
							y : yy(d.y),
							w : xx(d.z) - xx(0), 
							h : yy(d.y1) - yy(d.y),
							style : d.style, 
							data  : d.data, 
						};
					}))
					.enter().append('g')
						.attr("class", "data");
			} else if (type === 'heatplot') {
				var item = this.plotView.selectAll('.data')
					.data(mat.map(function(d) {return {
							x : xx(d.x),
							y : yy(d.y),
							w : xx(d.z) - xx(0), 
							h : yy(d.z) - yy(0), 
							style : d.style, 
							data  : d.data, 
						};
					}))
					.enter().append('g')
						.attr("class", "data");
			} else if (type === 'scatter') {
				var item = this.plotView.selectAll('.data')
					.data(mat.map(function(d) {return {
							x : xx(d.x),
							y : yy(d.y),
							r : 3*d.z, 
							style : d.style, 
							data  : d.data, 
						};
					}))
					.enter().append('g')
						.attr('class', 'data');
			} 
			plot_selected(item);
		} else {
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
		.call(d3.axisBottom(xx));
	this.plotView.append('g')
		.attr('transform', "translate(0, "+(self.plotDim.y0+self.plotDim.h)+")")
		.call(d3.axisTop(xx));
		
	this.plotView.append('g')
		.attr('transform', "translate("+self.plotDim.x0+",0)")
		.call(d3.axisLeft(yy));
	this.plotView.append('g')
		.attr('transform', "translate("(+self.plotDim.x0+self.plotDim.w)+",0)")
		.call(d3.axisRight(yy));

	return this;
};
